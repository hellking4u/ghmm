/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

/* This header can be read by both C and C++ compilers */
#ifndef _GHMM_XML_H
#define _GHMM_XML_H 1

#include <stdlib.h>
 
#ifdef __cplusplus

extern "C" {

typedef enum {
  DISCRETE,
  CONTINUOUS
} model_enum;

#include <ghmm/model.h>
struct model_wrapper {
  model_enum model_id;
  void * model_pt;
};
typedef struct model_wrapper model_t;

} /* extern C */

#else 

typedef enum {
  DISCRETE,
  CONTINUOUS
} model_enum;

#include <ghmm/model.h>
struct model_wrapper {
  model_enum model_id;
  void * model_pt;
};
typedef struct model_wrapper model_t;

#endif

/*
 * Prevent "name mangled" by C++ compiler
 * Declare a wrapper function to call C++ methods 
 */
#ifdef __cplusplus
extern "C" {
#endif
 
#if defined(__STDC__) || defined(__cplusplus)
  /* ANSI C prototypes */
  extern model_t* graphmldoc_cwrapper(char *filename);
#else
  /* K&R style */
  extern model_t* graphmldoc_cwrapper(char *filename);
#endif
 

#ifdef __cplusplus
}
#endif
 
#endif /*  end of header file */
