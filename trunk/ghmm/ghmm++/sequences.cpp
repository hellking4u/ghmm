/*
  author: Achim Gaedke
  created: 26. June 2001
  file: xmlio/examples/ghmm/sequence_reader.cpp
  $Id$
 */

#include <iostream>
#include <cerrno>
#include "sequences.h"
#include <xmlio/XMLIO_SkipObject.h>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

/***********************************************************************************/
/* templates are better! */

double_sequence::double_sequence(XMLIO_Attributes &attributes)
{
  XMLIO_Attributes::iterator id_key=attributes.find("id");
  if (id_key!=attributes.end())
    {
      id=id_key->second;
    }
  label=NULL;
}

double_sequence::~double_sequence()
{
  SAFE_DELETE(label);
}

XMLIO_Object* double_sequence::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag=="label")
    {
      if (label==NULL)
	{
	  label=new XMLIO_StringObject();
	  return label;
	}
      else
	{
	  cout<<"Only one label allowed in sequence"<<endl;
	}
    }
  else
    {
      cout<<tag<<" not supported in sequence"<<endl;
      return new XMLIO_SkipObject(this);
    }
}

/***********************************************************************************/

int_sequence::int_sequence(XMLIO_Attributes &attributes)
{
  XMLIO_Attributes::iterator id_key=attributes.find("id");
  if (id_key!=attributes.end())
    {
      id=id_key->second;
    }
  label=NULL;
}

int_sequence::~int_sequence()
{
  SAFE_DELETE(label);
}

XMLIO_Object* int_sequence::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag=="label")
    {
      if (label==NULL)
	{
	  label=new XMLIO_StringObject();
	  return label;
	}
      else
	{
	  cout<<"Only one label allowed in sequence"<<endl;
	}
    }
  else
    {
      cout<<tag<<" not supported in sequence"<<endl;
      return new XMLIO_SkipObject(this);
    }
}

int* int_sequence::create_array()
{
  int* array=(int*)malloc(sizeof(int)*size());
  XMLIO_ArrayObject<int>::iterator iter=begin();
  int i=0;
  while (iter!=end())
    {
      array[i]=*iter;
      i++;iter++;
    }
  return array;
}

int int_sequence::get_label()
{
  if (label!=NULL)
    {
      return strtol(label->c_str(),NULL,0);
    }
  else
    {
      return -1;
    }
      
}

double int_sequence::get_id()
{
  if (id.empty())
    return -1.0;
  else
    return strtod(id.c_str(),NULL);
}

/***********************************************************************************/

sequences_DiscretePD::sequences_DiscretePD(XMLIO_Attributes& attributes, const string& sequence_type)
{
  type=sequence_type;
  XMLIO_Attributes::iterator default_weight_key=attributes.find("default_weight");
  if (default_weight_key!=attributes.end())
    {
      errno=0;
      default_weight=strtod(default_weight_key->second.c_str(),NULL);
      if (errno)
	{
	  cout<<"default_weight is not a float!"<<endl;
	  default_weight=1;
	}
    }
  else
    default_weight=1;
  actual_weight=default_weight;
}

sequences_DiscretePD::~sequences_DiscretePD()
{
  while (!int_sequence_vector.empty())
    {
      SAFE_DELETE(double_sequence_vector.back());
      double_sequence_vector.pop_back();
    }
  while (!double_sequence_vector.empty())
    {
      SAFE_DELETE(double_sequence_vector.back());
      double_sequence_vector.pop_back();
    }
}

void sequences_DiscretePD::XMLIO_finishedReading()
{}

XMLIO_Object* sequences_DiscretePD::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag=="sequence")
    {
      if (type=="int")
	{
	  int_sequence* new_sequence=new int_sequence(attributes);
	  int_sequence_vector.push_back(new_sequence);
	  weight_vector.push_back(actual_weight);
	  actual_weight=default_weight;
	  return new_sequence;
	}
      else if (type=="double")
	{
	  double_sequence* new_sequence=new double_sequence(attributes);
	  double_sequence_vector.push_back(new_sequence);
	  weight_vector.push_back(actual_weight);
	  actual_weight=default_weight;
	  return new_sequence;
	}
      else
	{
	  cout<<"Do not know what to do with "<<tag<<" of type "<<type<<endl;
	  return new XMLIO_SkipObject(this);
	}      
    }
  else
    {
      cout<<"Do not know what to do with "<<tag<<endl;
      return new XMLIO_SkipObject(this);
    }
}

void sequences_DiscretePD::XMLIO_endTag(const string& tag)
{}

void sequences_DiscretePD::XMLIO_getCharacters(const string& characters)
{
  /* find floats */
  errno=0;
  actual_weight=strtod(characters.c_str(),NULL);
  if (errno)
    {
      cout<<"weight is not a float: "<<characters<<endl;
      actual_weight=default_weight;
    }
}

sequence_t* sequences_DiscretePD::create_sequence_t()
{
  sequence_t* seq=sequence_calloc(weight_vector.size());
  vector<double>::iterator weight_iter=weight_vector.begin();
  vector<int_sequence*>::iterator seq_iter=int_sequence_vector.begin();
  int i=0;

  while (weight_iter!=weight_vector.end() && seq_iter!=int_sequence_vector.end())
    {
      seq->seq_len[i]=(*seq_iter)->size();
      seq->seq[i]=(*seq_iter)->create_array();
      seq->seq_label[i]=(*seq_iter)->get_label();
      seq->seq_id[i]=(*seq_iter)->get_id();
      seq->seq_w[i]=*weight_iter;
      weight_iter++; seq_iter++; i++;
    }
  seq->seq_number=i;
  return seq;
}

sequence_d_t* sequences_DiscretePD::create_sequence_d_t()
{
  /* todo */
  return NULL;
}

/***********************************************************************************/

sequences::sequences()
{
  sequence_array=NULL;
}

sequences::sequences(XMLIO_Attributes& attributes)
{
  sequence_array=NULL;
  XMLIO_Attributes::iterator type_key=attributes.find("type");
  if (type_key!=attributes.end())
    {
      type=type_key->second;
      if (type_key->second!="int" && type_key->second!="double")
	{
	  cout<<"unsupported sequences type: "<<type_key->second<<endl;
	}
    }
  else
    cout<<"type not found"<<endl;
}

sequences::~sequences()
{
  SAFE_DELETE(sequence_array);
}

void sequences::XMLIO_finishedReading()
{
  /* prepare structures for return */
  cout<<"Done "<<toString()<<": "
      <<"Coding: "<<coding<<", "
      <<"Type: "<<type
      <<endl;
}

sequence_t* sequences::create_sequence_t() const
{
  sequence_t* int_seq=sequence_array->create_sequence_t();
  /* additional information add here */

  return int_seq;
}

sequence_d_t* sequences::create_sequence_d_t() const
{
  sequence_d_t* double_seq=sequence_array->create_sequence_d_t();
  /* additional information add here */
  return double_seq;
}

const char* sequences::toString() const
{
  return "sequences";
}

XMLIO_Object* sequences::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag=="DiscretePD")
    {
      if (sequence_array==NULL)
	{
	  sequence_array=new sequences_DiscretePD(attributes,type);
	  return sequence_array;
	}
      else
	{
	  cout<<"only one sequence array allowed in sequences"<<endl;
	}
    }
  else
    {
      cout<<tag<<" not allowed in sequences"<<endl;
      return new XMLIO_SkipObject(this);
    }
}

void sequences::XMLIO_endTag(const string& tag)
{}

/***********************************************************************************/

sequenceReader::sequenceReader()
{
  next_sequence_array=NULL;
}

sequenceReader::~sequenceReader()
{
  SAFE_DELETE(next_sequence_array);
  while (!empty())
    {
      SAFE_DELETE(front());
      pop_front();
    }

}

XMLIO_Object* sequenceReader::XMLIO_startTag(const string& tag, XMLIO_Attributes& attributes) {
  if (tag == "sequences") {
    next_sequence_array=new sequences(attributes);
    return next_sequence_array;
  }
  return new XMLIO_SkipObject(this);
}

void sequenceReader::XMLIO_endTag(const string& tag)
{
  if (tag=="sequences" && next_sequence_array!=NULL)
    {
      push_back(next_sequence_array);
      next_sequence_array=NULL;
    }
}

size_t sequenceReader::read_sequences(const string& filename)
{
  open(filename);
  read();
  close();
  return size();
}

/***********************************************************************************/

sequence_t* return_sequences(const char* file)
{
  sequenceReader* reader=new sequenceReader();
  (void)reader->read_sequences(file);
  sequence_t* seq=reader->front()->create_sequence_t();
  return seq;
}










