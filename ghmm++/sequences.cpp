/*
  author: Achim Gaedke
  created: 26. June 2001
  file: ghmm/ghmm++/sequences.cpp
  $Id$
 */

#include <cmath>
#include <iostream>
#include <cerrno>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "sequences.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

/***********************************************************************************/
/* templates are better! */

double_sequence::double_sequence(double* seq_data, size_t length)
{
  tag="sequence";
  for(size_t count=0; count<length; count++)
    push_back(seq_data[count]);
  label=NULL;
  id="";
}

double_sequence::double_sequence(int* seq_data, size_t length)
{
  tag="sequence";
  for(size_t count=0; count<length; count++)
    push_back(seq_data[count]);
  label=NULL;
  id="";
}

double_sequence::double_sequence(const string& name, XMLIO_Attributes &attributes)
{
  tag=name;
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

XMLIO_Element* double_sequence::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag=="label")
    {
      if (label==NULL)
	{
	  label=new XMLIO_StringElement();
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
    }
  return NULL;
}


const XMLIO_Attributes&
double_sequence::XMLIO_getAttributes () const{
  XMLIO_Attributes& attrs=(XMLIO_Attributes&)attributes;
  attrs.clear();
  if (!id.empty())
    attrs["id"]=id;
  return attributes;
}

const int double_sequence::XMLIO_writeContent (XMLIO_Document& doc) const
{
  int result=0;
  int tmp_result;
  if (label!=NULL) {
    tmp_result=doc.writeElement(*label);
    if (tmp_result<0) return tmp_result;
  }
  tmp_result=XMLIO_ArrayElement<double>::XMLIO_writeContent(doc);
  if (tmp_result<0) return tmp_result;
  return result+tmp_result;
}


double* double_sequence::create_double_array() const
{
  double* array=(double*)malloc(sizeof(double)*size());
  XMLIO_ArrayElement<double>::const_iterator iter=begin();
  int i=0;
  while (iter!=end())
    {
      array[i]=*iter;
      i++;iter++;
    }
  return array;
}

int* double_sequence::create_int_array() const
{
  int* array=(int*)malloc(sizeof(int)*size());
  XMLIO_ArrayElement<double>::const_iterator iter=begin();
  int i=0;
  while (iter!=end())
    {
      array[i]=(int)*iter;
      i++;iter++;
    }
  return array;
}

int double_sequence::get_label_as_int() const
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

double double_sequence::get_id_as_double() const
{
  if (id.empty())
    return -1.0;
  else
    return strtod(id.c_str(),NULL);
}

void double_sequence::print() const
{
  XMLIO_ArrayElement<double>::print();
}

/***********************************************************************************/

int_sequence::int_sequence(int* seq_data, size_t length)
{
  tag="sequence";
  for(size_t count=0; count<length; count++)
    push_back(seq_data[count]);
  label=NULL;
  id="";
}

int_sequence::int_sequence(double* seq_data, size_t length)
{
  tag="sequence";
  for(size_t count=0; count<length; count++)
    push_back((int)floor(seq_data[count]));
  label=NULL;
  id="";
}

int_sequence::int_sequence(const string& name, XMLIO_Attributes &attributes)
{
  tag=name;
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

XMLIO_Element* int_sequence::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag=="label")
    {
      if (label==NULL)
	{
	  label=new XMLIO_StringElement();
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
      return NULL;
    }
  return NULL;
}

const XMLIO_Attributes&
int_sequence::XMLIO_getAttributes () const{
  XMLIO_Attributes& attrs=(XMLIO_Attributes&)attributes;
  attrs.clear();
  if (!id.empty())
    attrs["id"]=id;
  return attributes;
}

const int int_sequence::XMLIO_writeContent (XMLIO_Document& doc) const
{
  int result=0;
  int tmp_result;
  if (label!=NULL) {
    tmp_result=doc.writeElement(*label);
    if (tmp_result<0) return tmp_result;
  }
  tmp_result=XMLIO_ArrayElement<int>::XMLIO_writeContent(doc);
  if (tmp_result<0) return tmp_result;
  return result+tmp_result;
}

int* int_sequence::create_int_array() const
{
  int* array=(int*)malloc(sizeof(int)*size());
  XMLIO_ArrayElement<int>::const_iterator iter=begin();
  int i=0;
  while (iter!=end())
    {
      array[i]=*iter;
      i++;iter++;
    }
  return array;
}

double* int_sequence::create_double_array() const
{
  double* array=(double*)malloc(sizeof(double)*size());
  XMLIO_ArrayElement<int>::const_iterator iter=begin();
  int i=0;
  while (iter!=end())
    {
      array[i]=(double)*iter;
      i++;iter++;
    }
  return array;
}


int int_sequence::get_label_as_int() const
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

double int_sequence::get_id_as_double() const
{
  if (id.empty())
    return -1.0;
  else
    return strtod(id.c_str(),NULL);
}

void int_sequence::print() const
{
  XMLIO_ArrayElement<int>::print();
}


/***********************************************************************************/


sequences_DiscretePD::sequences_DiscretePD(int** data, double* weight_data, size_t length, size_t number) {
  tag="DiscretePD";
  type="int";
  default_weight=1.0;
  for (size_t sequence_pos=0;sequence_pos<number;sequence_pos++)
    {
      int_sequence_vector.push_back(new int_sequence(data[sequence_pos],length));
      if (weight_data!=NULL)
	weight_vector.push_back(weight_data[sequence_pos]);
      else
	weight_vector.push_back(default_weight);
    }
}

sequences_DiscretePD::sequences_DiscretePD(double** data, double* weight_data, size_t length, size_t number) {
  tag="DiscretePD";
  type="double";
  default_weight=1.0;
  for (size_t sequence_pos=0;sequence_pos<number;sequence_pos++)
    {
      double_sequence_vector.push_back(new double_sequence(data[sequence_pos],length));
      if (weight_data!=NULL)
	weight_vector.push_back(weight_data[sequence_pos]);
      else
	weight_vector.push_back(default_weight);
    }
}

sequences_DiscretePD::sequences_DiscretePD(sequence_t* seq) {
  tag="DiscretePD";
  if (seq==NULL) return;
  type="int";
  default_weight=1.0;
  for (int sequence_pos=0;sequence_pos<seq->seq_number;sequence_pos++)
    {
      /* vector of doubles */
      int_sequence_vector.push_back(new int_sequence(seq->seq[sequence_pos],seq->seq_len[sequence_pos]));
      /* weight */
      weight_vector.push_back(seq->seq_w[sequence_pos]);
      /* label: missing todo */
      /* id: missing todo */
    }
}

sequences_DiscretePD::sequences_DiscretePD(sequence_d_t* seq) {
  tag="DiscretePD";
  if (seq==NULL) return;
  type="double";
  default_weight=1.0;
  for (int sequence_pos=0;sequence_pos<seq->seq_number;sequence_pos++)
    {
      /* vector of doubles */
      double_sequence_vector.push_back(new double_sequence(seq->seq[sequence_pos],seq->seq_len[sequence_pos]));
      /* weight */
      weight_vector.push_back(seq->seq_w[sequence_pos]);
      /* label: missing todo */
      /* id: missing todo */
    }
}


sequences_DiscretePD::sequences_DiscretePD (XMLIO_Attributes& attributes, const string& sequence_type) {
  tag="DiscretePD";
  type=sequence_type;
  XMLIO_Attributes::iterator default_weight_key=attributes.find("default_weight");
  if (default_weight_key!=attributes.end())
    {
      errno=0;
      default_weight=strtod(default_weight_key->second.c_str(),NULL);
      if (errno)
	{
	  cout<<"default_weight is not a float!"<<endl;
	  default_weight=1.0;
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
      SAFE_DELETE(int_sequence_vector.back());
      int_sequence_vector.pop_back();
    }
  while (!double_sequence_vector.empty())
    {
      SAFE_DELETE(double_sequence_vector.back());
      double_sequence_vector.pop_back();
    }
}

void sequences_DiscretePD::XMLIO_finishedReading()
{}

XMLIO_Element* sequences_DiscretePD::XMLIO_startTag(const string& name, XMLIO_Attributes &attributes)
{
  if (name=="sequence")
    {
      if (type=="int")
	{
	  int_sequence* new_sequence=new int_sequence(name,attributes);
	  int_sequence_vector.push_back(new_sequence);
	  weight_vector.push_back(actual_weight);
	  actual_weight=default_weight;
	  return new_sequence;
	}
      else if (type=="double")
	{
	  double_sequence* new_sequence=new double_sequence(name,attributes);
	  double_sequence_vector.push_back(new_sequence);
	  weight_vector.push_back(actual_weight);
	  actual_weight=default_weight;
	  return new_sequence;
	}
      else
	{
	  cout<<"Do not know what to do with "<<name<<" of type "<<type<<endl;
	  return NULL;
	}      
    }
  else
    {
      cout<<"Do not know what to do with "<<name<<endl;
      return NULL;
    }
}

void sequences_DiscretePD::XMLIO_endTag(const string& tag)
{}

void sequences_DiscretePD::XMLIO_getCharacters(const string& characters)
{
  /* find weight floats */
  errno=0;
  const char* old_pos=characters.c_str();
  while (*old_pos!=0 && isspace(*old_pos)) old_pos++;
  if (*old_pos==0) 
    {
      actual_weight=default_weight;
      return;
    }

  char* new_pos=(char*)old_pos;
  actual_weight=strtod(old_pos,&new_pos);
  if (errno || new_pos==old_pos)
    {
      cout<<"weight is not a float: "<<characters<<endl;
      actual_weight=default_weight;
    }
}

const int sequences_DiscretePD::XMLIO_writeContent (XMLIO_Document& doc) const {
  int result = doc.writeEndl();
  if (result<0) return result;
  /* distinguish between int and double sequences */
  if (type=="int") {
      vector<double>::const_iterator weight_pos=weight_vector.begin();
      for (vector<int_sequence*>::const_iterator pos=int_sequence_vector.begin();
	   pos != int_sequence_vector.end();
	   ++pos)
	{
	  int tmp_result;
	  if (weight_pos != weight_vector.end())
	    {
	      tmp_result=doc.writef("%f ",*weight_pos);
	      if (tmp_result<0) return tmp_result;
	      result+=tmp_result;
	      ++weight_pos;
	    }
	  tmp_result=doc.writeElement(**pos);
	  if (tmp_result<0) return tmp_result;
	  result+=tmp_result+1;
	  tmp_result=doc.writeEndl();
	  if (tmp_result<0) return tmp_result;
	}
      return result;
  }
  if (type=="double") {
      vector<double>::const_iterator weight_pos=weight_vector.begin();
      for (vector<double_sequence*>::const_iterator pos=double_sequence_vector.begin();
	   pos != double_sequence_vector.end();
	   ++pos)
	{
	  int tmp_result;
	  if (weight_pos != weight_vector.end())
	    {
	      tmp_result=doc.writef("%f ",*weight_pos);
	      if (tmp_result<0) return tmp_result;
	      result+=tmp_result;
	      ++weight_pos;
	    }
	  tmp_result=doc.writeElement(**pos);
	  if (tmp_result<0) return tmp_result;
	  result+=tmp_result+1;
	  tmp_result=doc.writeEndl();
	  if (tmp_result<0) return tmp_result;
	}
      return result;
  }
}


const XMLIO_Attributes& sequences_DiscretePD::XMLIO_getAttributes () const {
  XMLIO_Attributes& attrs=(XMLIO_Attributes&)attributes;
  attrs.clear();
  if (default_weight!=0) {
    strstream my_stream;
    my_stream<<default_weight<<ends;
    attrs["default_weight"]=my_stream.str();
  }
  return attributes;
}

sequence_t* sequences_DiscretePD::create_sequence_t() const
{
  sequence_t* seq=sequence_calloc(weight_vector.size());
  vector<double>::const_iterator weight_iter=weight_vector.begin();
  vector<int_sequence*>::const_iterator seq_iter=int_sequence_vector.begin();
  int i=0;

  while (weight_iter!=weight_vector.end() && seq_iter!=int_sequence_vector.end())
    {
      seq->seq_len[i]=(*seq_iter)->size();
      seq->seq[i]=(*seq_iter)->create_int_array();
      seq->seq_label[i]=(*seq_iter)->get_label_as_int();
      seq->seq_id[i]=(*seq_iter)->get_id_as_double();
      seq->seq_w[i]=*weight_iter;
      weight_iter++; seq_iter++; i++;
    }
  seq->seq_number=i;
  return seq;
}

sequence_d_t* sequences_DiscretePD::create_sequence_d_t() const
{
  sequence_d_t* seq=sequence_d_calloc(weight_vector.size());
  vector<double>::const_iterator weight_iter=weight_vector.begin();
  vector<double_sequence*>::const_iterator seq_iter=double_sequence_vector.begin();
  int i=0;

  while (weight_iter!=weight_vector.end() && seq_iter!=double_sequence_vector.end())
    {
      seq->seq_len[i]=(*seq_iter)->size();
      seq->seq[i]=(*seq_iter)->create_double_array();
      seq->seq_label[i]=(*seq_iter)->get_label_as_int();
      seq->seq_id[i]=(*seq_iter)->get_id_as_double();
      seq->seq_w[i]=*weight_iter;
      weight_iter++; seq_iter++; i++;
    }
  seq->seq_number=i;
  return seq;
}

void sequences_DiscretePD::print() const
{
  if (type=="int")
    {
      cout<<"type is "<<type<<endl;
      vector<double>::const_iterator weight_pos=weight_vector.begin();
      for (vector<int_sequence*>::const_iterator pos=int_sequence_vector.begin();
	   pos != int_sequence_vector.end();
	   ++pos)
	{
	  if (weight_pos != weight_vector.end())
	    {
	      cout<<*weight_pos<<" ";
	      ++weight_pos;
	    }
	  (*pos)->print();
	}
    }
  else if (type=="double")
    {
      cout<<"type is "<<type<<endl;
      vector<double>::const_iterator weight_pos=weight_vector.begin();
      for (vector<double_sequence*>::const_iterator pos=double_sequence_vector.begin();
	   pos != double_sequence_vector.end();
	   ++pos)
	{
	  if (weight_pos != weight_vector.end())
	    {
	      cout<<*weight_pos<<" ";
	      ++weight_pos;
	    }
	  (*pos)->print();
	}      
    }
  else
    {
      cerr<<toString()<<": type "<<type<<" not supported!"<<endl;
    }
}

/***********************************************************************************/

sequences::sequences()
{
  tag="sequences";
  type="";
  sequence_array=NULL;
}

sequences::sequences(int** data, double* weights, size_t length, size_t number)
{
  tag="sequences";
  type="int";
  sequence_array=new sequences_DiscretePD(data,weights,length,number); 
}

sequences::sequences(double** data, double* weights, size_t length, size_t number)
{
  tag="sequences";
  type="double";
  sequence_array=new sequences_DiscretePD(data,weights,length,number);
}

sequences::sequences(sequence_t* seq)
{
  tag="sequences";
  type="int";
  sequence_array=new sequences_DiscretePD(seq);
}

sequences::sequences(sequence_d_t* seq)
{
  tag="sequences";
  type="double";
  sequence_array=new sequences_DiscretePD(seq);  
}

sequences::sequences(const string& name, XMLIO_Attributes& attributes)
{
  tag=name;
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
  return;
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

void sequences::print() const
{
  cout<<toString()<<":"<<endl;
  if (sequence_array!=NULL)
    sequence_array->print();
  else
    cout<<"no sequences available"<<endl;
}

XMLIO_Element* sequences::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
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
    }
  return NULL;
}

void sequences::XMLIO_endTag(const string& tag)
{}

const XMLIO_Attributes&
sequences::XMLIO_getAttributes () const{
  XMLIO_Attributes& attrs=(XMLIO_Attributes&)attributes;
  attrs.clear();
  if (!type.empty())
    attrs["type"]=type;
  if (!coding.empty())
    attrs["coding"]=coding;
  return attributes;
}

bool
sequences::XMLIO_isEmpty () const{
  return false;
}

const int
sequences::XMLIO_writeContent (XMLIO_Document& doc) const{
  if (sequence_array!=NULL) {
    doc.writeEndl();
    return doc.writeElement(*sequence_array);
    doc.writeEndl();
  }
  else {
    return doc.writeComment("no sequences");
  }
}








