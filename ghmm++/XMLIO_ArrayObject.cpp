#include "XMLIO_ArrayObject.h"

/***********************************************************************************/

void XMLIO_StringObject::XMLIO_getCharacters(const string& characters)
{
  append(characters);
}

/***********************************************************************************/

int* get_token(const string &characters,int &position)
{
  position++;
  return new int(0);
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
