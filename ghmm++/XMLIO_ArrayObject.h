#ifndef _XMLIO_ARRAYOBJECT_H
#define _XMLIO_ARRAYOBJECT_H 1

#include <iostream>
#include <vector>
#include <string>
#include <xmlio/XMLIO_Object.h>

int* get_token(const string&, int&);

template<class T>
class XMLIO_Array: public XMLIO_Object, public vector<T> {
 public:
  void print() const{
    XMLIO_Array<T>::const_iterator pos=begin();
    while (pos!=end()) {cout<<*pos<<endl;pos++;}
  }

  void XMLIO_getCharacters(const string& characters){
    T* new_member;
    int position=0;
    do
      {
	new_member=get_token(characters, position);
	if (new_member!=NULL)
	  {
	    push_back(*new_member);
	    SAFE_DELETE(new_member);
	  }
      } 
    while (position<characters.size());
  }

  char word_seperator;
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
