#include "XMLIO_ArrayObject.h"

#include <string>
#include <cstdlib>
#include <cerrno>
#include <cstdio>

/***********************************************************************************/

void XMLIO_StringObject::XMLIO_getCharacters(const string& characters)
{
  append(characters);
}

/***********************************************************************************/

size_t XMLIO_evaluate_token(const string &characters, const size_t position, const size_t length, int* &return_int)
{
  char* new_pos;
  const char* old_pos=&(characters.c_str()[position]);
  errno=0;
  int new_value=strtol(old_pos,&new_pos,0);
  if (errno || old_pos==new_pos)
    {
      return_int=NULL;
      return 0;
    }
  else
    {
      return_int=new int(new_value);
      return new_pos-old_pos;
    }
}

size_t XMLIO_evaluate_token(const string &characters, const size_t position, const size_t length, double* &return_double)
{
  char* new_pos;
  const char* old_pos=&(characters.c_str()[position]);
  errno=0;
  double new_value=strtod(old_pos,&new_pos);
  if (errno || old_pos==new_pos)
    {
      return_double=NULL;
      return 0;
    }
  else
    {
      return_double=new double(new_value);
      return new_pos-old_pos;
    }
}


size_t XMLIO_evaluate_token(const string &characters, const size_t position, const size_t length, string* &return_string)
{
  return_string=new string(characters,position,length);
  return length;
}

size_t XMLIO_evaluate_token(const string &characters, const size_t position, const size_t length, char* &return_char)
{
  return_char=new char(characters[position]);
  return 1;
}
