#ifndef _XMLIO_ARRAYOBJECT_H
#define _XMLIO_ARRAYOBJECT_H 1

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <xmlio/XMLIO_Object.h>

size_t evaluate_token(const string& characters , const size_t position, const size_t length,int* &return_int);
size_t evaluate_token(const string& characters , const size_t position, const size_t length, double* &return_double);

template<class T>
class XMLIO_Array: public XMLIO_Object, public vector<T> {
 public:
  void print() const{
    XMLIO_Array<T>::const_iterator pos=begin();
    while (pos!=end()) {cout<<*pos<<endl;pos++;}
  }

  void XMLIO_getCharacters(const string& characters){
    for (int position=0 ; position<characters.size(); position++)
      {
	/* delete leading spaces */
	if (isspace(characters[position])) continue;
	int length=0;
	while (position+length<characters.size() && !isspace(characters[length])) length++;
	T* new_element;
	(void)evaluate_token(characters, position, length, new_element);
	if (new_element!=NULL)
	  push_back(*new_element);
	else
	  cerr<<"can not evaluate "<<characters.substr(position,length)<<endl;
	SAFE_DELETE(new_element);
	position+=length;
      }
  }
};

class XMLIO_StringObject: public XMLIO_Object, public string
{
 public:
  void XMLIO_getCharacters(const string& characters);
};


class XMLIO_IntArrayObject: public XMLIO_Object, public vector<int>
{
 public:
  void XMLIO_getCharacters(const string& characters);
};

class XMLIO_DoubleArrayObject: public XMLIO_Object, public vector<double>
{
 public:
  void XMLIO_getCharacters(const string& characters);
};

#endif /* _XMLIO_ARRAYOBJECT_H */
