/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/ghmm.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ghmm.h"
#include <xmlio/XMLIO_IgnoreObject.h>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

ghmm::ghmm(const string& tag, XMLIO_Attributes &attributes)
{
  /* do some initialisation */ 
  hmm_model=NULL;
  shmm_model=NULL;
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
  return_model->s=(state*)malloc(sizeof(state)*return_model->N); /* state pointer array */
  if (return_model->s==NULL)
    {
      cerr<<"could not allocate states array"<<endl;
      return NULL;
    }

  /* initialise states */
  int state_counter=0;
  /* iterates over existing nodes */
  vector<gxl_node*>::const_iterator node_pos=my_graph->vector<gxl_node*>::begin();
  while (node_pos!=my_graph->vector<gxl_node*>::end())
    {
      const string& node_id=(*node_pos)->id;
      state* this_state=&(return_model->s[state_counter]);
      /* emission probabilities */
      this_state->b=(double*)malloc(sizeof(double)*return_model->M);
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


      size_t array_pos; /* position in double array of state */
      set<int>::const_iterator tranisiton_pos; /* position in set of transitons */

      /* transitions to_from */
      /* the inverse adiascence list */
      const set<int>& to_from_transition_idx=my_graph->get_to_from_transitions(node_id);
      this_state->in_states=to_from_transition_idx.size(); /* number of incoming states */
      this_state->in_id=(int*)malloc(sizeof(int)*this_state->in_states); /* id of incoming states */
      this_state->in_a=(double*)malloc(sizeof(double)*this_state->in_states); /* prob of incoming states */
      array_pos=0;
      tranisiton_pos=to_from_transition_idx.begin();
      while (tranisiton_pos!=to_from_transition_idx.end())
	{
	  this_state->in_id[array_pos]=*tranisiton_pos;
	  ghmm_edge* my_edge=dynamic_cast<ghmm_edge*>(((vector<gxl_edge*>)*my_graph)[*tranisiton_pos]);
	  if (my_edge!=NULL && !my_edge->empty())
	    {
	      this_state->in_a[array_pos]=my_edge->front();
	    }
	  else
	    this_state->in_a[array_pos]=0;
	  ++tranisiton_pos;
	  array_pos++;
	}

      /* transitions from_to */
      /* adiascense list */
      const set<int>& from_to_transition_idx=my_graph->get_from_to_transitions(node_id);
      this_state->out_states=from_to_transition_idx.size(); /* number of incoming states */
      this_state->out_id=(int*)malloc(sizeof(int)*this_state->out_states); /* id of incoming states */
      this_state->out_a=(double*)malloc(sizeof(double)*this_state->out_states); /* prob of incoming states */
      array_pos=0;
      tranisiton_pos=from_to_transition_idx.begin();
      while (tranisiton_pos!=from_to_transition_idx.end())
	{
	  this_state->out_id[array_pos]=*tranisiton_pos;
	  ghmm_edge* my_edge=dynamic_cast<ghmm_edge*>(((vector<gxl_edge*>)*my_graph)[*tranisiton_pos]);
	  if (my_edge!=NULL && !my_edge->empty())
	    {
	      this_state->out_a[array_pos]=my_edge->front();
	    }
	  else
	    this_state->out_a[array_pos]=0;
	  ++tranisiton_pos;
	  array_pos++;
	}

      /* emission probabilities */
      this_state->b=(double*)calloc(sizeof(double),return_model->M);
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
	    }
	} /* ghmm_Emissions!=NULL */

      ++node_pos;
      state_counter++;
    }

  return return_model;
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



