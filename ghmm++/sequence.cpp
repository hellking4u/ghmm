/*
  author: Achim Gaedke
  created: 2001-09-08
  file: ghmm/ghmm++/sequence.cpp
  $Id$
 */

#include <cmath>
#include <iostream>
#include <cerrno>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ghmm++/sequence.h"
#include "ghmm/sequence.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

/***********************************************************************************/

double_sequence::double_sequence(const string& name, XMLIO_Attributes &attrs)
  :sequence<double>(name, attrs){
}

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

/***********************************************************************************/

int_sequence::int_sequence(const string& name, XMLIO_Attributes &attrs)
  :sequence<int>(name, attrs){
}


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

