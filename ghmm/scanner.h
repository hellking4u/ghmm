/*********************************************************************
  author       : Frank Nübel
  filename     : ghmm/ghmm/scanner.h
  created      : TIME: 11:42:57     DATE: Tue 27. August 1996
  $Id$

__copyright__

*********************************************************************/

#ifndef SCANNER_H
#define SCANNER_H

#include <stdio.h>

#ifdef __cplusplus 
extern "C" {
#endif

/**
   @name scanner structure and methods
*/

/*@{ */

/*
 */
typedef struct scanner_t {
  FILE* fp;              /* Parameter-File, das eingelesen wird        */
  char* filename;        /* Name dieses Files                          */
  int   line;            /* current line                               */
  int   pos;             /* position in the current text line          */
  int   idlen;           /* maximal length of the identifier           */
  char* id;              /* identifier                                 */
  int   txtlen;          /* maximal length of the text line            */
  char* txt;             /* contains the current line text  
			    (used for error message )                  */
  char  c;               /* current char                               */
  char  esc;             /* current char is escaped by '\'             */
  char  err;             /* err == 0 : OK                              */
  char  eof;             /* eof == 1 : end of file reached             */
  
  /*******************/
  int   resolution_used; /* is set after the first use of length units */
  float x_resolution;    /* x-size of one dot in inch;                 */
  float y_resolution;    /* y-size of one dot in inch;                 */
  float x_scale;
  float y_scale;
} scanner_t;

  /**
   */
scanner_t* scanner_alloc( char* filename );
  /**
   */
int        scanner_consume( scanner_t* s, char ch );
  /**
   */
int        scanner_consume_block( scanner_t* s );
  /**
   */
int        scanner_error( scanner_t* s, char* message );
  /**
   */
int        scanner_free( scanner_t** s );
  /**
   */
int        scanner_free_array( int* len, void***arr );

  /**
   */
void*      scanner_get_array( scanner_t* s, int*len, char* type );
  /**
   */
double     scanner_get_double( scanner_t* s );
  /**
   */
int        scanner_get_id( scanner_t* s );
  /**
   */
int        scanner_get_int( scanner_t* s );
  /**
   */
int        scanner_get_name( scanner_t* s );
  /**
   */
char*      scanner_get_path( scanner_t* s );
  /**
   */
char*      scanner_get_str( scanner_t* s, int* len, int cmode );
  /**
   */
double**   scanner_get_d_matrix(scanner_t* s, int*rows, int*cols);

/**************************/
  /**
   */
int        scanner_get_index(scanner_t* s, int n);
  /**
   */
int        scanner_get_length_x(scanner_t* s);
  /**
   */
int        scanner_get_length_y(scanner_t* s);
  /**
   */
double     scanner_get_resolution(scanner_t* s);

  /**
   */
#define scanner_get_boolean( s )         (!!scanner_get_int( s ))
  /**
   */
#define scanner_get_char( s )            ((char)(scanner_get_int( s )))
  /**
   */
#define scanner_get_cchar( s )           ((char)(scanner_get_int( s )))
  /**
   */
#define scanner_get_cstring( s )         scanner_get_str( (s), NULL, 1 )
  /**
   */
#define scanner_get_string( s )          scanner_get_str( (s), NULL, 0 )

  /**
   */
#define scanner_get_char_array(s,len)    scanner_get_array((s),(len),"char" )
  /**
   */
#define scanner_get_cstring_array(s,len) scanner_get_array((s),(len),"cstring" )
  /**
   */
#define scanner_get_double_array(s,len)  scanner_get_array((s),(len),"double" )
  /**
   */
#define scanner_get_double_earray(s,len) scanner_get_array((s),(len),"edouble" )
  /**
   */
#define scanner_get_int_array(s,len)     scanner_get_array((s),(len),"int" )
  /**
   */
#define scanner_get_string_array(s,len)  scanner_get_array((s),(len),"string" )


  /*@} scanner section */

#ifdef __cplusplus 
}
#endif

#endif /* SCANNER_H */
