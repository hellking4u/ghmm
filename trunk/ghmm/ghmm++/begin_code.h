/*
  created: 31 Jan 2002 by Peter Pipenbacher
  authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
  file   : $Source$
  $Id$

  __copyright__
  
*/


/* This shouldn't be nested -- included it around code only. */
#ifdef _begin_code_h
#error Nested inclusion of begin_code.h
#endif
#define _begin_code_h

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

/* Windows has namespaces */
#if defined(WIN32) 
#  ifdef _MSC_VER
#    pragma warning(disable: 4251)
#  endif
#  if ! defined(HAVE_NAMESPACES)
#    define HAVE_NAMESPACES
#  endif
#endif

/* Some compilers use a special export keyword */
#ifndef DECLSPEC
# ifdef __BEOS__
#  if defined(__GNUC__)
#   define DECLSPEC	__declspec(dllexport)
#  else
#   define DECLSPEC	__declspec(export)
#  endif
# else
# ifdef WIN32
#  ifdef GHMM_EXPORTS
#    define DECLSPEC	__declspec(dllexport)
#  else
#    define DECLSPEC	__declspec(dllimport)
#  endif
# else
#  define DECLSPEC
# endif
# endif
#endif

/* some standard includes */
#ifdef HAVE_CSTRING
#  include <cstring>
#else
#  include <string.h>
#endif
#ifdef HAVE_CSTDLIB
#  include <cstdlib>
#else
#  include <stdlib.h>
#endif
#ifdef HAVE_CSTDIO
#  include <cstdio>
#else
#  include <stdio.h>
#endif
#ifdef HAVE_CMATH
#  include <cmath>
#else
#  include <math.h>
#endif

/* define MIN, MAX macros */
#define MAX(a,b) ((a > b) ? a : b)
#define MIN(a,b) ((a < b) ? a : b)

#if defined(HAVE_NAMESPACES) && defined(__SUNPRO_CC)
/* for SC5.0 in /usr/include/floatingpoint.h (SunOs 2.7) */
using std::FILE;
#endif
