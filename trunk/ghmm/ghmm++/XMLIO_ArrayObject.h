#ifndef _XMLIO_ARRAYOBJECT_H
#define _XMLIO_ARRAYOBJECT_H 1

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <xmlio/XMLIO_Definitions.h>
#include <xmlio/XMLIO_Object.h>
#include <xmlio/XMLIO_SkipObject.h>

size_t XMLIO_evaluate_token(const string& characters , const size_t position, const size_t length,int* &return_int);
size_t XMLIO_evaluate_token(const string& characters , const size_t position, const size_t length, double* &return_double);
size_t XMLIO_evaluate_token(const string& characters , const size_t position, const size_t length, string* &return_string);
size_t XMLIO_evaluate_token(const string& characters , const size_t position, const size_t length, char* &return_char);

template<class T>
class XMLIO_ArrayObject: public XMLIO_Object, public vector<T> {
 public:
  void print() const{
    XMLIO_ArrayObject<T>::const_iterator pos=begin();
    while (pos!=end()) {cout<<*pos<<endl;pos++;}
  }

  void XMLIO_getCharacters(const string& characters){
    for (int position=0 ; position<characters.size(); position++)
      {
	/* delete leading spaces */
	if (isspace(characters[position])) continue;
	/* look forward to end of token */
	int new_position=position+1;
	while (new_position<characters.size() && !isspace(characters[new_position]))
	  new_position++;
	/* evaluate token */
	T* new_element;
	(void)XMLIO_evaluate_token(characters, position, new_position-position, new_element);
	if (new_element!=NULL)
	  {
	    /* copy it to vector */
	    push_back(*new_element);
	    SAFE_DELETE(new_element);
	  }
	else
	  {
	    cerr<<"can not evaluate "<<characters.substr(position, new_position-position)
		<<" between "<<position<<" and "<<new_position<<endl;

	  }
	/* jump behind it and the next white space (or beyond end) */
	position=new_position;
      }
  }
};

class XMLIO_StringObject: public XMLIO_Object, public string
{
 public:
  void XMLIO_getCharacters(const string& characters);
};

template<class E>
class XMLIO_ElementArray: public vector<E>, public XMLIO_Object
{
  /* to do */
};

template<class C, class E>
class XMLIO_ContentElementPair
{
 public:
  C content;
  E element;
};

template<class C, class E>
class XMLIO_ContentElementPairObject: public vector<XMLIO_ContentElementPair<C,E> >, public XMLIO_Object
{
 public:
  void XMLIO_getCharacters(const string& characters)
    {
      return;
    }
  XMLIO_Object* XMLIO_startElement()
    {
      return new XMLIO_SkipObject(this); 
    }
  void XMLIO_endElement()
    {}
};


#endif /* _XMLIO_ARRAYOBJECT_H */
