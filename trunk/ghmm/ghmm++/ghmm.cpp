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

hmm::hmm(const string& tag, XMLIO_Attributes &attributes)
{
  /* do some initialisation */ 
  hmm_model=NULL;
  shmm_model=NULL;
  Initial=NULL;
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

}

XMLIO_Object* hmm::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  /* what's next? */
  if (tag=="Graph")
    {
      cerr<<toString()<<": "<<tag<<" not implemented... DoItNOW"<<endl;
      return new XMLIO_IgnoreObject();
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
  else if (tag=="Emissions")
    {
      cerr<<toString()<<": "<<tag<<" not implemented... DoItNOW"<<endl;
      return new XMLIO_IgnoreObject();      
    }
  else
    {
      cerr<<toString()<<": found unexpected element "<<tag<<", ignoring"<<endl;
      return new XMLIO_IgnoreObject();
    }
}

void hmm::XMLIO_endTag(const string& tag)
{
  /* not needed now */
}

void hmm::XMLIO_getCharacters(const string& characters)
{
  /* do nothing...*/
}

void hmm::XMLIO_finishedReading()
{
  /* sanity check */
  cerr<<"sanity check ToDoNOW"<<endl;
}


void hmm::print() const
{
  cout<<toString()<<": "<<"ToDoNOW!"<<endl;
  if (Initial!=NULL)
      Initial->print();
  else
      cout<<"InitialStates empty"<<endl;
  cout<<"Graph: ToDoNOW"<<endl
      <<"Emissions: ToDoNOW"<<endl
      <<"Alphabet: ToDoNOW"<<endl;
}

const char* hmm::toString() const
{
  return "hmm";
}
