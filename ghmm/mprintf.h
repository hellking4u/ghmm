/*******************************************************************************
  author       : Frank Nübel
  filename     : ghmm/ghmm/mprintf.h
  created      : TIME: 11:27:32     DATE: Wed 14. May 1997
  $Id$

__copyright__

*******************************************************************************/

#ifndef MPRINTF_H
#define MPRINTF_H

#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

///  
char* mprintf( char* dst, int maxlen, char* format, ... );
///  
char* mprintf_dyn( char* dst, int maxlen, char* format, ... );
///  
char* mprintf_va( char* dst, int maxlen, char* format, va_list args );
///  
char* mprintf_va_dyn( char* dst, int maxlen, char* format, va_list args );



#ifdef __cplusplus
}
#endif /*__cplusplus*/
#endif /* MPRINTF_H */
