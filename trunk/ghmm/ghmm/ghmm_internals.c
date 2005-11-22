/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/logging.c
*       Authors:  Janne Grunau
*
*       Copyright (C) 1998-2004 Alexander Schliep 
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*	Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik, 
*                               Berlin
*                                   
*       Contact: schliep@ghmm.org             
*
*       This library is free software; you can redistribute it and/or
*       modify it under the terms of the GNU Library General Public
*       License as published by the Free Software Foundation; either
*       version 2 of the License, or (at your option) any later version.
*
*       This library is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*       Library General Public License for more details.
*
*       You should have received a copy of the GNU Library General Public
*       License along with this library; if not, write to the Free
*       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*
*       This file is version $Revision$ 
*                       from $Date$
*             last change by $Author$.
*
*******************************************************************************/


#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ghmm.h"
#include "ghmm_internals.h"


static char * qmessage;

static int maxlevel = 3;

static void (* logfunc)(int level, const char * message, void * clientdata);

static void * logfunc_data;


void ighmm_logging (int level, const char * proc, const char * str) {

  char * message;
  int len = strlen(proc);

  /* concatenate the precompiler info with the actual logging text if any */
  if (str || qmessage) {
    
    /* get the queued messege text */
    if (!str) {
      str = qmessage;
      qmessage = NULL;
    }

    len += strlen(str);
    
    message = malloc(sizeof(char)*(len+1));
    if (!message)
      return;
    
    message = strcpy(message, proc);
    message = strcat(message, str);
    
  }
  else
    message = (char *)proc;

  /* if defined use external logging function */
  if (logfunc) {
    logfunc(level, message, logfunc_data);
  }
  /* otherwise simmple logging stderr */
  else 
    if (level < maxlevel) {
      switch (level) {
      case LDEBUG:
	fputs("DEBUG: ", stderr);
	break;
      case LINFO:
	fputs("INFO: ", stderr);
	break;
      case LWARN:
	fputs("WARNING: ", stderr);
	break;
      case LERROR:
	fputs("ERROR: ", stderr);
	break;
      case LCRITIC:
	fputs("CRITICAL: ", stderr);
	break;
      default:
	break;
      }
      fputs(message, stderr);
    }
}

void ighmm_queue_mes(char * text) {

  /* let no message get lost */
  if (!qmessage)
    qmessage = text;
  else {
    ighmm_logging(LCRITIC, __FILE__":ighmm_queue_mes line "TOSTRING(__LINE__)": unable to queue message since it exists an unprocessed message! ", text);
    exit(3);
  }
}

void ghmm_set_logfunc(void (* fptr)(int, const char *, void *), void * clientdata) {

  logfunc = fptr;
  logfunc_data = clientdata;

}

void ghmm_set_loglevel(int level) {

  maxlevel = level;

}
