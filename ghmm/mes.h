/*********************************************************************
  author       : Frank Nübel
  filename     : /vol/bspk/src/hmm/lib/mes.h
  created      : TIME: 17:33:09     DATE: Thu 25. January 1996
  last-modified: TIME: 15:05:26     DATE: Thu 29. January 1998
*********************************************************************/

#ifndef MES_H
#define MES_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

#define MES_0_PTR 0
#define MES_NEG_ARG 1
#define MES_BIG_ARG 2
#define MES_0_ARG 3
#define MES_FALSE_ARG 4

#define MES_FLAG_TIME        1
#define MES_FLAG_WIN         4
#define MES_FLAG_FILE       16 
#define MES_FLAG_FILE_WIN   ( MES_FLAG_FILE | MES_FLAG_WIN )
#define MES_FLAG_TIME_WIN   ( MES_FLAG_TIME | MES_FLAG_FILE | MES_FLAG_WIN )

#define MES_MODUL_INFO  "(" __DATE__ ":" __FILE__ ")"
#define MES_PROC_INFO   "(" __DATE__ ":" __FILE__  ":" CUR_PROC ")"


#define MES_WIN             MES_FLAG_WIN, -1, NULL, NULL
#define MES_FILE            MES_FLAG_FILE, -1, NULL, NULL
#define MES_TIME            (MES_FLAG_TIME | MES_FLAG_FILE), -1, NULL, NULL
#define MES_FILE_WIN        MES_FLAG_FILE_WIN, -1, NULL, NULL
#define MES_TIME_WIN        MES_FLAG_TIME_WIN, -1, NULL, NULL
#define MES_PROT_WIN        MES_FLAG_WIN, __LINE__ ,MES_PROC_INFO, CUR_PROC
#define MES_PROT_FILE       MES_FLAG_FILE, __LINE__ ,MES_PROC_INFO, CUR_PROC
#define MES_PROT            MES_FLAG_FILE_WIN, __LINE__ ,MES_PROC_INFO, CUR_PROC
#define MES_PROT_TIME       MES_FLAG_TIME_WIN, __LINE__ ,MES_PROC_INFO, CUR_PROC

///  
#define mes_check_ptr( arg, call) \
  if(!(arg)) { mes_err( #arg, MES_0_PTR, MES_PROC_INFO ); call; } else
///  
#define mes_check_bnd( arg, bnd, call ) \
  if((arg) > (bnd)){ mes_err( #arg, MES_BIG_ARG, MES_PROC_INFO ); call;} else
///  
#define mes_check_0( arg, call ) \
  if(!(arg)) { mes_err( #arg, MES_0_ARG, MES_PROC_INFO ); call;} else
///  
#define mes_check_neg( arg, call ) \
  if((arg)<0){ mes_err( #arg, MES_NEG_ARG, MES_PROC_INFO );call;} else
///  
#define mes_check_expr( arg, call ) \
  if(!(arg)){ mes_err( #arg, MES_FALSE_ARG, MES_PROC_INFO );call;} else
  
  
///  
#define mes_file( txt )     mes_smart( MES_FLAG_FILE, txt, -1 )
///  
#define mes_file_win( txt ) mes_smart( MES_FLAG_FILE_WIN, txt, -1 )
///  
#define mes_win( txt )      mes_smart( MES_FLAG_WIN, txt, -1 )
///  
#define mes_clock_abs( txt1, txt2 ) mes( MES_WIN, "%s(%T)%s", txt1, txt2 )
///  
#define mes_clock_rel( txt1, txt2 ) mes( MES_WIN, "%s(%t)%s", txt1, txt2 )
  
///  
#define mes_proc()      mes( MES_PROT, NULL )
///  
#define mes_prot( txt ) mes( MES_PROT_TIME, txt )

///
char* mes_get_std_path( void );

///
void  mes( int flags, int line, char* xproc, char* proc, char* format, ... );
///
void  mes_smart( int flags, char* txt, int bytes );
///
int   mes_ability( int on );
///
void  mes_err( char* txt, int error_nr, char* proc_info );
///
void  mes_exit( void );
///
void  mes_fformat( char* txt, char* filename, int line, char* proc_info );
///
void  mes_init( char* logfile, void(*winfct)(char*), int argc, char* argv[] );
///
void  mes_init_args( int argc, char*argv[]);
///
void  mes_init_logfile( char*file_name );
///
void  mes_init_winfct( void(win_fct)(char*) );
///
int   mes_insert( FILE* fp, char src, int cnt );
///
void  mes_time( void );
///
void  mes_va( int flags,int line,char*xproc,char*proc,char*format,va_list args);
///
int   mes_win_ability( int on );

///
void* mes_calloc( int bytes );
///
int   mes_copy( char *oldname, char *newname );
///
int   mes_fgetc( FILE*fp );
///
int   mes_fflush( FILE* fp );
///
FILE* mes_fopen( char* filename, char* attribute_string );
///
int   mes_fprintf( FILE* fp, char* format, ... );
///
int   mes_fputc( FILE*fp, char chr );
///
int   mes_fputs( FILE*fp, char* str );
///
int   mes_fread( FILE* fp, void* mem, int bytes );
///
int   mes_fread_quiet( FILE* fp, void* mem, int bytes );
///
int   mes_fseek( FILE*fp, long offset, int fromwhere );
#ifdef WIN32
///
int mes_fseek64( FILE *fp, unsigned int uoff, unsigned int loff, int fromwhere );
#endif
///
int   mes_ftell( FILE* fp );
///
int   mes_fwrite( FILE* fp, void* mem, int bytes );
///
void* mes_malloc( int bytes );
///
int   mes_move( char* oldname, char* newname );
///
int   mes_realloc( void** mem, int bytes );
///
int   mes_remove( char* filename );
///
int   mes_rename( char* oldname, char* newname );
///
FILE* mes_tmpfile( void );
#if 0
///
char* mes_tmpname(char *str);
#endif

#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif /* MES_H */
