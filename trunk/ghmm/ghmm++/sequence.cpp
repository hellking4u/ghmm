/*
  author: Achim Gaedke
  created: 2001-09-08
  file: ghmm/ghmm++/sequence.cpp
  $Id$
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <strstream>
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
  label="";
  id="";
}

double_sequence::double_sequence(int* seq_data, size_t length)
{
  tag="sequence";
  for(size_t count=0; count<length; count++)
    push_back(seq_data[count]);
  label="";
  id="";
}

double_sequence::double_sequence(sequence_d_t* seq, int sequence_pos) {
  tag="sequence";
  int length = seq->seq_len[sequence_pos];
  double* seq_data=seq->seq[sequence_pos];
  for(int count=0; count < length; count++)
    push_back(seq_data[count]);
  strstream tmpstr;
  strstream tmpstr2;
  tmpstr<<seq->seq_label[sequence_pos]<<ends;
  label = tmpstr.str(); tmpstr.clear();
  tmpstr2.setf(ios::fixed);
  tmpstr2<<setprecision(0)<<seq->seq_id[sequence_pos]<<ends;
  id = tmpstr2.str();
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


/***********************************************************************************/

int_sequence::int_sequence(const string& name, XMLIO_Attributes &attrs)
  :sequence<int>(name, attrs){
}


int_sequence::int_sequence(int* seq_data, size_t length)
{
  tag="sequence";
  for(size_t count=0; count<length; count++)
    push_back(seq_data[count]);
  label="";
  id="";
}

int_sequence::int_sequence(double* seq_data, size_t length)
{
  tag="sequence";
  for(size_t count=0; count<length; count++)
    push_back((int)floor(seq_data[count]));
  label="";
  id="";
}

int_sequence::int_sequence(sequence_t* seq, int sequence_pos) {
  tag="sequence";
  int length = seq->seq_len[sequence_pos];
  int* seq_data=seq->seq[sequence_pos];
  for(int count=0; count < length; count++)
    push_back(seq_data[count]);
  strstream tmpstr;
  strstream tmpstr2;
  tmpstr<<seq->seq_label[sequence_pos]<<ends;
  label = tmpstr.str(); tmpstr.clear();
  tmpstr2.setf(ios::fixed);
  tmpstr2<<setprecision(0)<<seq->seq_id[sequence_pos]<<ends;
  id = tmpstr2.str();
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


