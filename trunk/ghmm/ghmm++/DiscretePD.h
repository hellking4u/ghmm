
#include <iostream>
#include <vector>
#include <xmlio/XMLIO_Object.h>
#include <xmlio/XMLIO_ObjectReader.h>
#include <xmlio/XMLIO_SkipObject.h>
#include <xmlio/XMLIO_IgnoreObject.h>
#include <xmlio/XMLIO_ArrayObject.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  template <class E>
    class DiscretePD: public XMLIO_ContentElementArrayObject<double,E>
    {
    public:
      /**
	 if no weight is given default_weight is 1.
       */
      DiscretePD()
	{
	  default_weight=1;
	  actual_weight=default_weight;
	}

      /**
	 looks for attribute default_weight
       */
      DiscretePD(const string& tag, XMLIO_Attributes &attributes)
	{
	  XMLIO_Attributes::iterator default_weight_key=attributes.find("default_weight");
	  if (default_weight_key!=attributes.end())
	    {
	      double* tmp_default;
	      (void)XMLIO_evaluate_token(pos->second,0,pos->second->size(),&tmp_default);
	      if (tmp_default==NULL)
		{
		  cerr<<toString()<<": default weight is not a floating point value"<<endl;
		  default_weight=1;
		}
	      else
		{
		  default_weight=*tmp_default;
		  delete tmp_default;
		}
	    }
	  else
	    default_weight=1;
	  actual_weight=default_weight;
	}
    private:
      double default_weight;
      double actual_weight;
    };

  class State: public XMLIO_Object
    {
    public:
      State(const string& tag, XMLIO_Attributes &attributes);
      void print() const;
      const char* toString() const;
    private:
      string id_ref;
    };
  
  class InitialStates: public XMLIO_Object
    {
    public:
      InitialStates(const string& tag, XMLIO_Attributes &attributes);
      XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      const char* toString() const {return "InitialStates";}
      void print() const;
    private:
      DiscretePD<State>* state_pd;      
    };

  /**
     Parses all elements of a xml file, optional filters elements by name and appends them to a vector
   */

  template <class E>
  class XMLIO_ElementArrayReader: public XMLIO_ObjectReader, public vector<E*>
    {
    public:
      XMLIO_ElementArrayReader(): vector<E*>()
	{
	  element_name="";
	}
      
      ~XMLIO_ElementArrayReader()
	{
	  while (!empty())
	    {
	      SAFE_DELETE(back());
	      pop_back();
	    }
	}

      XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
	{
	  cerr<<"found element "<<tag<<endl;
	  if (element_name.empty() || element_name==tag)
	    {
	      E* actual_element=new E(tag, attributes);
	      push_back(actual_element);
	      return actual_element;
	    }
	  else
	    return new XMLIO_SkipObject(this);
	}

      /**
	 This is the name of the parsed and saved elements
       */
      void set_element_name(const string& new_name)
	{element_name=new_name;}

      const char* toString()
	{return "XMLIO_ElementArrayReader";}

      /**
	 Reads the elements from a xml file
       */
      size_t read_file(const char* filename)
	{
	  open(filename);
	  read();
	  close();
	  return size();
	}

    private:
      /**
	 Saves the name of the elements to parse
       */
      string element_name;
    };

#ifdef HAVE_NAMESPACES
}
#endif



