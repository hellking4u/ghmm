/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/ghmm.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ghmm/rng.h>
#include <ghmm/vector.h>
#include <ghmm/sequence.h>
#include "ghmm++/ghmm.h"
#include <xmlio/XMLIO_IgnoreObject.h>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

ghmm::ghmm(const string& tag, XMLIO_Attributes &attributes)
{
  /* do some initialisation */ 
  Initial=NULL;
  ghmm_Emissions=NULL;
  my_graph=NULL;
  type='\0';

  /* read attributes */
  
  XMLIO_Attributes::const_iterator pos;
  /* madatory argument */  
  pos=attributes.find("type");
  if (pos!=attributes.end())
    {
      if (pos->second=="discrete")
	{
	  /* discrete type */
	  type='d';
	}
      else if (pos->second=="continuous")
	{
	  /* continuous type */
	  type='c';
	}
      else if (pos->second=="switched")
	{
	  /* switched type */
	  type='s';
	}
      else
	{
	  cerr<<toString()<<": unknown type "<<pos->second<<endl;
	}
    }
  else
    {
      /**/
      cerr<<toString()<<": type attribute missing"<<endl;
    }

  pos=attributes.find("prior");
  if (pos!=attributes.end())
    {
      double* tmp_prior_pointer=NULL;
      (void)XMLIO_evaluate_token(pos->second,0,pos->second.size(),tmp_prior_pointer);
      if (tmp_prior_pointer==NULL)
	{
	  cerr<<toString()<<": attribute value of prior is not a double"<<endl;
	  prior=1;
	}
      else
	{
	  prior=*tmp_prior_pointer;
	  delete tmp_prior_pointer;
	}
    }
}

XMLIO_Object* ghmm::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  /* what's next? */
  if (tag=="graph")
    {
      if (my_graph!=NULL)
	{
	  cerr<<toString()<<" only one graph expected"<<endl;
	  return NULL;
	}
      else
	{
	  my_graph=new ghmm_graph(tag,attributes);
	  return my_graph;
	}
    }
  else if (tag=="InitialStates")
    {
      if (Initial!=NULL)
	{
	  cerr<<toString()<<": initial States allready existing, appending new section"<<endl;
	}
      else
	{
	  Initial=new InitialStates(tag,attributes);
	}
      return Initial;
    }
  else if (tag=="Alphabet")
    {
      if (my_alphabet==NULL)
	{
	  my_alphabet=new ghmm_alphabet(tag,attributes);
	  return my_alphabet;
	}
      else
	{
	  cerr<<toString()<<": alphabet allready existing,ignoring"<<endl;
	  return new XMLIO_IgnoreObject();
	}

    }
  else if (tag=="Emissions")
    {
      if (ghmm_Emissions!=NULL)
	{
	  cerr<<toString()<<" only one Emissions section expected"<<endl;
	  return new XMLIO_IgnoreObject();
	}
      else
	{
	  ghmm_Emissions=new Emissions(tag,attributes);
	  return ghmm_Emissions;
	}
    }
  else
    {
      cerr<<toString()<<": found unexpected element "<<tag<<", ignoring"<<endl;
      return new XMLIO_IgnoreObject();
    }
}

void ghmm::XMLIO_endTag(const string& tag)
{
  /* not needed now */
}

void ghmm::XMLIO_getCharacters(const string& characters)
{
  /* do nothing...*/
}

void ghmm::XMLIO_finishedReading()
{
  /* done at graph...*/
  /* more needed here */
}

void ghmm::print() const
{
  cout<<toString()<<":"<<endl;
  if (Initial!=NULL)
      Initial->print();
  else
      cout<<"InitialStates empty"<<endl;
  if (my_graph!=NULL)
    my_graph->print();
  else
    cout<<"no graph!"<<endl;
  if (ghmm_Emissions!=NULL)
    ghmm_Emissions->print();
  else
    cout<<"no emissions!"<<endl;
  if (my_alphabet!=NULL)
      my_alphabet->print();
  else
      cout<<"alphabet not specified"<<endl;
}

const char* ghmm::toString() const
{
  return "ghmm";
}

model* ghmm::create_model() const
{
  /* model type mismatch */
  if (type!='d')
    return (model*)NULL;

  /* are necessary informations available? */
  if (my_graph==NULL || Initial==NULL)
    {
      return (model*)NULL;
    }

  /* alpahabet and emissions */

  /* allocate model */
  model* return_model=(model*)malloc(sizeof(model));
  if (return_model==NULL)
    {
      cerr<<"could not allocate model structure"<<endl;
      return (model*)NULL;
    }
  
  return_model->N=my_graph->vector<gxl_node*>::size();   /* number of states */
  return_model->prior=prior;

  /* number of symbols */
  if (my_alphabet!=NULL)
    {
      return_model->M=my_alphabet->size();
    }
  else
    {
      return_model->M=0;
    }

  /* allocate state array */
  return_model->s=(state*)calloc(return_model->N,sizeof(state)); /* state pointer array */
  if (return_model->s==NULL)
    {
      cerr<<"could not allocate states array"<<endl;
      return NULL;
    }

  /* initialise states */
  int state_counter=0;
  /* iterates over existing nodes */
  vector<gxl_node*>::const_iterator node_pos=my_graph->vector<gxl_node*>::begin();
  /* for initial state possibilty normalization */
  double initial_state_prob_sum=0;
  while (node_pos!=my_graph->vector<gxl_node*>::end())
    {
      const string& node_id=(*node_pos)->id;
      state* this_state=&(return_model->s[state_counter]);
      /* emission probabilities */
      this_state->b=(double*)calloc(return_model->M,sizeof(double));
      if (this_state->b==NULL)
	{
	  cerr<<"could not allocate emmision probabilities vector"<<endl;
	  return NULL;
	}

      /* search for Initial State probability */
      map<State*,double>* initial_state_map=Initial->get_map();
      map<State*,double>::iterator state_pos=initial_state_map->begin();
      while(state_pos!=initial_state_map->end() && (state_pos->first)->get_id()!=node_id)
	++state_pos;
      if (state_pos!=initial_state_map->end())
	this_state->pi=state_pos->second;
      else
	this_state->pi=0;
      /* for normalization */
      initial_state_prob_sum+=this_state->pi;

      size_t array_pos; /* position in double array of state */
      set<int>::const_iterator tranisiton_pos; /* position in set of transitons */

      /* transitions from_to */
      /* adjacent list */
      const set<int>& from_to_transition_idx=my_graph->get_from_to_transitions(node_id);
      this_state->out_states=from_to_transition_idx.size(); /* number of incoming states */
      this_state->out_id=(int*)calloc(this_state->out_states,sizeof(int)); /* id of incoming states */
      this_state->out_a=(double*)calloc(this_state->out_states,sizeof(double)); /* prob of incoming states */
      array_pos=0;
      tranisiton_pos=from_to_transition_idx.begin();
      while (tranisiton_pos!=from_to_transition_idx.end())
	{
	  ghmm_edge* my_edge=dynamic_cast<ghmm_edge*>(((vector<gxl_edge*>)*my_graph)[*tranisiton_pos]);
	  if (my_edge!=NULL)
	    {
	      /* look for to state index */
	      this_state->out_id[array_pos]=my_graph->get_node_idx(my_edge->to);
	      /* look for weight */
	      if (!my_edge->empty())
		{
		  this_state->out_a[array_pos]=my_edge->front();
		}
	      else
		this_state->out_a[array_pos]=0;
	    }
	  else
	    {
	      cerr<<toString()<<": From State "<<state_counter<<": can't convert to ghmm Edge: Dynamic cast failed!"<<endl;
	    }
	  ++tranisiton_pos;
	  array_pos++;
	}
      /* normalize out_states */
      vector_normalize(this_state->out_a,this_state->out_states);

      /* emission probabilities */
      this_state->b=(double*)calloc(return_model->M,sizeof(double));
      if (ghmm_Emissions!=NULL)
	{
	  /* find emission information */
	  Emissions::const_iterator emissions_pos=ghmm_Emissions->begin();
	  while(emissions_pos!=ghmm_Emissions->end() && (*emissions_pos)->state!=node_id)
	    ++emissions_pos;
	  if (emissions_pos!=ghmm_Emissions->end())
	    {
	      const Emission& my_emission=**emissions_pos;
	      /* found it */
	      /* set fix state */
	      int emission_prob_idx=0;
	      this_state->fix=my_emission.fix;
	      /* iterate over sequence of probabilities */
	      for (Emission::const_iterator emission_pos=my_emission.begin();
		   emission_pos!=my_emission.end() && emission_prob_idx<return_model->M;
		   ++emission_pos)
		{
		  if (emission_pos->content!=NULL)
		      this_state->b[emission_prob_idx++]=*(emission_pos->content);
		} /* for emission_pos */
	      /* normalize emission probs */
	      vector_normalize(this_state->b,return_model->M);
	    }
	} /* ghmm_Emissions!=NULL */

      ++node_pos;
      state_counter++;
    }

  /* normalize initial State Distribution */
  if (initial_state_prob_sum!=0)
    for (int state_idx=0;state_idx<return_model->N;state_idx++)
      return_model->s[state_idx].pi/=initial_state_prob_sum;


  /* transitions to_from */
  /* the inverse adjacent list, derived from the adjacent list */

  /* collect information into map*/
  map<int, map<int, double> > inverse_adjacent;
  for(int state_idx=0; state_idx<return_model->N; state_idx++)
    {
      state& this_state=return_model->s[state_idx];
      for(int transition=0; transition<this_state.out_states; transition++)
	inverse_adjacent[this_state.out_id[transition]][state_idx]=this_state.out_a[transition];
    }

  /* put them into the structures */
  for(map<int, map<int, double> >::iterator state_pos=inverse_adjacent.begin(); state_pos!=inverse_adjacent.end(); ++state_pos)
    {
      state this_state=return_model->s[state_pos->first];
      /* allocate arrays */
      this_state.in_states=state_pos->second.size();
      this_state.in_id=(int*)calloc(this_state.in_states,sizeof(int));
      this_state.in_a=(double*)calloc(this_state.in_states,sizeof(double));
      /* write infos from map */
      int array_count=0;
      for(map<int, double>::iterator transition_pos=state_pos->second.begin(); transition_pos!=state_pos->second.end(); ++transition_pos)
	{
	 this_state.in_id[array_count]=transition_pos->first;
	 this_state.in_a[array_count]=transition_pos->second;
	 array_count++;
	}
    }

  return return_model;
}

sequences* ghmm::generate_sequences(int number, int length)
{
  /* create the model struct */
  model* my_model=create_model();
  if (my_model==NULL)
    return NULL;

  model_print(stdout,my_model);
  /* is it initialised before ? */
  gsl_rng_init();
  /* create some sequences */
  sequence_t* my_sequences=model_generate_sequences(my_model,0,length,number);
  /* and delete the model */
  model_free(&my_model);
  if (my_sequences==NULL)
      return NULL;

  /* make ghmm++ objects */
  sequences* return_seq=new sequences(my_sequences);

  sequence_free(&my_sequences);
  return return_seq;
}


const string& ghmm::get_id() const
{
  return id;
}

smodel* ghmm::create_smodel() const
{
  if (type!='c' || type!='s') return (smodel*)NULL;

  return NULL;
}


ghmm::~ghmm()
{
  SAFE_DELETE(ghmm_Emissions);
  SAFE_DELETE(my_graph);
  SAFE_DELETE(Initial);
}



