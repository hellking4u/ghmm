#include "XMLIO_ArrayObject.h"

#include <cstdlib>
#include <cerrno>
#include <cstdio>

/***********************************************************************************/

void XMLIO_StringObject::XMLIO_getCharacters(const string& characters)
{
  append(characters);
}

/***********************************************************************************/

size_t evaluate_token(const string &characters, const size_t position, const size_t length, int* &return_int)
{
  errno=0;
  int new_value=strtol(&(characters.c_str()[position]),NULL,0);
  if (errno)
    {
      return_int=NULL;
      return 0;
    }
  else
    {
      return_int=new int(new_value);
      return 1;
    }
}

size_t evaluate_token(const string &characters, const size_t position, const size_t length, double* &return_double)
{
  errno=0;
  double new_value=strtod(&(characters.c_str()[position]),NULL);
  if (errno)
    {
      return_double=NULL;
      return 1;
    }
  else
    {
      return_double=new double(new_value);
      return 0;
    }
}


/***********************************************************************************/

/* relies on strings, that are not truncated */
void XMLIO_IntArrayObject::XMLIO_getCharacters(const string& characters)
{

  int value;
  const char* start=characters.c_str();
  const char* pos=start;
  const char* end=&start[characters.size()];

#if 0
  fprintf(stderr,"found sequence data");
  fwrite(s,len,1,stderr);
#endif

  /* parse the data and save them */

  while (pos<end)
    {
      /* jump over it */
      if (isspace(*pos)) pos++;
      /* try to convert to integer */
      else
        {
          char* tmp;
          errno = 0;
          /* Parse it. */
          value = strtol (pos, &tmp, 0);
          /* Add it in, if not overflow.  */
          if (errno)
            fprintf(stderr,"value parsing error\n");
          else
            {
              /* store value */	      
	      push_back(value);
            }
	  /* next please! */
          pos=(const char*)tmp;
        }
    }
}

/***********************************************************************************/

/* relies on strings, that are not truncated */
void XMLIO_DoubleArrayObject::XMLIO_getCharacters(const string& characters)
{
  double value;
  const char* start=characters.c_str();
  const char* pos=start;
  const char* end=&start[characters.size()];

#if 0
  fprintf(stderr,"found sequence data");
  fwrite(s,len,1,stderr);
#endif

  /* parse the data and save them */

  while (pos<end)
    {
      /* jump over it */
      if (isspace(*pos)) pos++;
      /* try to convert to integer */
      else
        {
          char* tmp;
          errno = 0;
          /* Parse it. */
          value = strtod (pos, &tmp);
          /* Add it in, if not overflow.  */
          if (errno)
            fprintf(stderr,"value parsing error\n");
          else
            {
              /* store value */	      
	      push_back(value);
            }
	  /* next please! */
          pos=(const char*)tmp;
        }
    }
}
