#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "DiscretePD.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

State::State(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag=="State")
    {
      XMLIO_Attributes::const_iterator id_pos=attributes.find("id");
      if (id_pos==attributes.end())
	cerr<<"no id in State's attributes found!"<<endl;
      else
	id_ref=id_pos->second;
    }
  else
    {
      cerr<<tag<<": this is not State, that was expected."<<endl;
    }
}

void State::print() const
{
  cout<<"State "<<id_ref<<endl;
}

const char* State::toString() const
{
  return "State";
}

InitialStates::InitialStates(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag!="InitialStates")
    {
      cerr<<"content not expected"<<endl;
    }
  state_pd=NULL;  
}

XMLIO_Object* InitialStates::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag=="DiscretePD")
    {
      if (state_pd!=NULL)
	{
	  cerr<<"Only one DiscretePD section allowed"<<endl;
	  return new XMLIO_IgnoreObject();
	}
      else
	{
	  state_pd=new DiscretePD<State>;
	  /* only State elements are welcome */
	  state_pd->set_element_name("State");
	  return state_pd;
	}
    }
}

void InitialStates::print() const
{
  if (state_pd==NULL)
    {
      cout<<"no InitialState distribution"<<endl;
    }
  else
    {
      state_pd->print();
    }
}
