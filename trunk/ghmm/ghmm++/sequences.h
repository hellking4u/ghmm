#ifndef __GHMMPP_SEQUENCES__
#define __GHMMPP_SEQUENCES__

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <xmlio/XMLIO_Definitions.h>
#include <xmlio/XMLIO_Object.h>
#include <xmlio/XMLIO_StringObject.h>
#include <xmlio/XMLIO_ObjectReader.h>
#include <xmlio/XMLIO_ArrayObject.h>
#include <ghmm/sequence.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

sequence_t* std::return_sequences(const char* file);

class double_sequence: public XMLIO_ArrayObject<double>
{
public:
  double_sequence(XMLIO_Attributes &attributes);
  ~double_sequence();
  XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
  double* create_double_array() const;
  int* create_int_array() const;
  int get_label_as_int() const;
  double get_id_as_double() const;
private:
  XMLIO_StringObject* label;
  string id;
};

class int_sequence: public XMLIO_ArrayObject<int>
{
public:
  int_sequence(XMLIO_Attributes &attributes);
  ~int_sequence();
  XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
  double* create_double_array() const;
  int* create_int_array() const;
  int get_label_as_int() const;
  double get_id_as_double() const;
private:
  XMLIO_StringObject* label;
  string id;
};

class sequences_DiscretePD: public XMLIO_Object
{
public:
  sequences_DiscretePD(XMLIO_Attributes& attributes, const string& type);
  ~sequences_DiscretePD();
  XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
  void XMLIO_endTag(const string& tag);
  void XMLIO_finishedReading();
  void XMLIO_getCharacters(const string& characters);
  const string& get_type() const {return type;}

  sequence_t* create_sequence_t() const;
  sequence_d_t* create_sequence_d_t() const;

private:
  string type;
  double actual_weight, default_weight;
  vector<int_sequence*> int_sequence_vector;
  vector<double_sequence*> double_sequence_vector;
  vector<double> weight_vector;
};

class sequences: public XMLIO_Object
{
public:
  /**/
  sequences();
  sequences(XMLIO_Attributes& attributes);
  ~sequences();
  XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
  void XMLIO_endTag(const string& tag);
  void XMLIO_finishedReading();
  const char* toString() const;
  const string& get_type() const {return type;}
  sequence_t* create_sequence_t() const;
  sequence_d_t* create_sequence_d_t() const;

private:
  sequences_DiscretePD* sequence_array;
  string type;
  string coding;
};

class sequenceReader: public XMLIO_ObjectReader, public list<sequences*>
{

public:
  sequenceReader();
  ~sequenceReader();
  size_t read_sequences(const string& filename);
  XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes& attributes);
  void XMLIO_endTag(const string& tag);

private:
  sequences* next_sequence_array;
};

#ifdef HAVE_NAMESPACES
}
#endif


#endif /* __GHMMPP_SEQUENCES__ */
