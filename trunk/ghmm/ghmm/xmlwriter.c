/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/xmlwriter.h
*       Authors:  Janne Grunau
*
*       Copyright (C) 1998-2006 Alexander Schliep 
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*	Copyright (C) 2002-2006 Max-Planck-Institut fuer Molekulare Genetik, 
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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

#include "ghmm.h"
#include "mes.h"
#include "mprintf.h"
#include "ghmm_internals.h"
#include "xmlwriter.h"


#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

#define MY_ENCODING "ISO-8859-1"


static xmlChar *ConvertInput(const char *in, const char *encoding);


/* ========================================================================= */
static int writeIntAttribute(xmlTextWriterPtr writer, const char * name, int value) {
#define CUR_PROC "writeIntAttribute"
  int rc;
  char * estr, * str = ighmm_mprintf(NULL, 0, "%d", value);
  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST name,
				   BAD_CAST str);
  m_free(str);
  if (rc < 0) {
    estr = ighmm_mprintf(NULL, 0, "Failed to write attribute %s (%d)", name, value);
    ighmm_queue_mes(estr);
    return -1;
  }
  else return 0;
#undef CUR_PROC
}

/* ========================================================================= */
static int writeDoubleAttribute(xmlTextWriterPtr writer, const char * name, double value) {
#define CUR_PROC "writeDoubleAttribute"
  int rc;
  char * estr, * str = ighmm_mprintf(NULL, 0, "%.8g", value);
  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST name,
				   BAD_CAST str);
  m_free(str);
  if (rc < 0) {
    estr = ighmm_mprintf(NULL, 0, "Failed to write attribute %s (%d)", name, value);
    ighmm_queue_mes(estr);
    return -1;
  }
  else return 0;
#undef CUR_PROC
}

/* ========================================================================= */
static char * strModeltype(int modelType) {
#define CUR_PROC "strModelType"

  char * mt;

  ARRAY_CALLOC(mt, 200);
  
  if (modelType > 0) {
    if (modelType & GHMM_kLeftRight)
      strcat(mt, "left-right ");
    if (modelType & GHMM_kSilentStates)
      strcat(mt, "silent ");
    if (modelType & GHMM_kTiedEmissions)
      strcat(mt, "tied ");
    if (modelType & GHMM_kHigherOrderEmissions)
      strcat(mt, "higher-order ");
    if (modelType & GHMM_kBackgroundDistributions)
      strcat(mt, "background ");
    if (modelType & GHMM_kLabeledStates)
      strcat(mt, "labeled ");
    if (modelType & GHMM_kTransitionClasses)
      strcat(mt, "transition-classes ");
    if (modelType & GHMM_kDiscreteHMM)
      strcat(mt, "discrete ");
    if (modelType & GHMM_kContinuousHMM)
      strcat(mt, "continuous ");
    if (modelType & GHMM_kPairHMM)
      strcat(mt, "pair ");
  } else {
    GHMM_LOG(LERROR, "can't write models with unspecified modeltype");
    goto STOP;
  }
  
  return mt;
 STOP:
  m_free(mt);
  return NULL;
#undef CUR_PROC
}

/* ========================================================================= */
static char * doubleArrayToCSV(double * array, int size) {
#define CUR_PROC "doubleArrayToCSV"

  int i, pos=0;
  char * csv;
  int maxlength = (10+2)*size;

  ARRAY_MALLOC(csv, maxlength);

  for (i=0; i<size-1 && pos<maxlength-10; i++) {
    pos += sprintf(csv+pos, "%.8g, ", array[i]);
  }
  if (i<size-1) {
    GHMM_LOG(LERROR, "writing CSV failed");
    goto STOP;
  } else {
    pos += sprintf(csv+pos, "%.8g", array[i]);
  }
  /*printf("%d bytes of %d written\n", pos, maxlength);*/
  return csv;
STOP:
  return NULL;
#undef  CUR_PROC
}


/* ========================================================================= */
static int writeAlphabet(xmlTextWriterPtr writer, alphabet_s * alfa) {
#define CUR_PROC "writeAlphabet"

  int i, rc;

  rc = xmlTextWriterStartElement(writer, BAD_CAST "alphabet");
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
    goto STOP;;
  }

  if (writeIntAttribute(writer, "id", alfa->id))
    GHMM_LOG_QUEUED(LERROR);

  for (i=0; i<alfa->size; i++) {
    rc = xmlTextWriterStartElement(writer, BAD_CAST "symbol");
    if (rc < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
      goto STOP;;
    }
    if (writeIntAttribute(writer, "code", i))
      GHMM_LOG_QUEUED(LERROR);

    rc = xmlTextWriterWriteRaw(writer, (xmlChar *)(alfa->symbols[i]));
    if (rc < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterWriteRaw");
      goto STOP;
    }

    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement");
      goto STOP;
    }
  }

  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement");
    goto STOP;
  }

  return 0;
STOP:
  return -1;
#undef CUR_PROC
}

/* ========================================================================= */
static int writeBackground(xmlTextWriterPtr writer, ghmm_d_background_distributions * bg) {
#define CUR_PROC "writeBackground"

  int i, rc;
  char * tmp;

  for (i=0; i<bg->n; i++) {
    rc = xmlTextWriterStartElement(writer, BAD_CAST "background");
    if (rc < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
      goto STOP;;
    }
    
    printf("background %d name: %s\n", i, bg->name[i]);
    if (xmlTextWriterWriteAttribute(writer, BAD_CAST "key", BAD_CAST (bg->name[i])))
      GHMM_LOG(LERROR, "Error at xmlTextWriterWriteAttribute");
      
    if (writeIntAttribute(writer, "order", bg->order[i]))
      GHMM_LOG_QUEUED(LERROR);
    
    tmp = doubleArrayToCSV(bg->b[i], pow(bg->m, bg->order[i]+1));
    if (tmp) {
      rc = xmlTextWriterWriteRaw(writer, BAD_CAST tmp);
      m_free(tmp);
      if (rc < 0) {
	GHMM_LOG(LERROR, "Error at xmlTextWriterWriteRaw");
	goto STOP;
      }
    } else {
      GHMM_LOG(LERROR, "converting array to CSV failed");
      m_free(tmp);
      goto STOP;
    }

    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement");
      goto STOP;
    }
  }
  
  return 0;
STOP:
  return -1;
#undef CUR_PROC
}

/* ========================================================================= */
static int writeDiscreteStateContents(xmlTextWriterPtr writer, fileData_s * f,
				      int moNo, int sNo) {
#define CUR_PROC "writeDiscreteStateContents"

  int bgId, cLabel, rc, order, tied;
  char * tmp=NULL;

  /* writing discrete distribution */
  rc = xmlTextWriterStartElement(writer, BAD_CAST "discrete");
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
    goto STOP;
  }

  if (f->model.d[moNo]->s[sNo].fix)
    if (writeIntAttribute(writer, "fixed", 1)) {
      GHMM_LOG_QUEUED(LERROR);
      goto STOP;
    }

  if ((f->model.d[moNo]->model_type & GHMM_kHigherOrderEmissions)
      && f->model.d[moNo]->order[sNo]) {
    order = f->model.d[moNo]->order[sNo];
    if (writeIntAttribute(writer, "order", order)) {
      GHMM_LOG_QUEUED(LERROR);
      goto STOP;
    }
  } else
    order = 0;

  tmp = doubleArrayToCSV(f->model.d[moNo]->s[sNo].b, pow(f->model.d[moNo]->M, order+1));
  if (tmp) {
    rc = xmlTextWriterWriteRaw(writer, BAD_CAST tmp);
    m_free(tmp);
    if (rc < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterWriteRaw");
      goto STOP;
    }
  } else {
    GHMM_LOG(LERROR, "converting array to CSV failed");
    m_free(tmp);
    goto STOP;
  }

  /* end discrete distribution */
  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement");
    goto STOP;
  }

  /* writing backgroung key */
  if (f->model.d[moNo]->model_type & GHMM_kBackgroundDistributions) {
    bgId = f->model.d[moNo]->background_id[sNo];
    if (bgId > -1) {
      if (f->model.d[moNo]->bp->name[bgId]) {
	rc = xmlTextWriterWriteElement(writer, BAD_CAST "backgroundKey",
				       BAD_CAST f->model.d[moNo]->bp->name[bgId]);
	if (rc<0) {
	  GHMM_LOG(LERROR, "Error at xmlTextWriterWriteElement");
	  goto STOP;
	}
      } else {
	GHMM_LOG(LERROR, "background name is NULL pointer, invalid model");
	goto STOP;
      }
    }
  }

  /* writing class label */
  if (f->model.d[moNo]->model_type & GHMM_kLabeledStates) {
    cLabel = f->model.d[moNo]->label[sNo];
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "class","%d", cLabel);
    if (rc<0) {
      GHMM_LOG(LERROR, "failed to write class label");
      goto STOP;
    }
  }
  
  /* duration (not implemented yet, maybe never */
#if 0
  if (f->model.d[moNo]->model_type & GHMM_kBackgroundDistributions) {
    bgId = f->model.d[moNo]->background_id[sNo];
    if (bgId > -1) {
      if (f->model.d[moNo]->bp->name[bgId]) {
	rc = xmlTextWriterWriteElement(writer, BAD_CAST "backgroundKey",
				       BAD_CAST f->model.d[moNo]->bp->name[bgId]);
	if (rc<0) {
	  GHMM_LOG(LERROR, "Error at xmlTextWriterWriteElement");
	  goto STOP;
	}
      } else {
	GHMM_LOG(LERROR, "background name is NULL pointer, invalid model");
	goto STOP;
      }
    }
  }
#endif

  /* writing positions */
  if ((f->model.d[moNo]->s[sNo].xPosition > 0)
      && (f->model.d[moNo]->s[sNo].xPosition > 0)) {

    if (xmlTextWriterStartElement(writer, BAD_CAST "position") < 0) {
      GHMM_LOG(LERROR, "failed to start position element"); goto STOP;}
    if (xmlTextWriterWriteFormatAttribute(writer, "x", "%d", f->model.d[moNo]->s[sNo].xPosition) < 0) {
      GHMM_LOG(LERROR, "failed to write x position"); goto STOP;    }
    if (xmlTextWriterWriteFormatAttribute(writer, "y", "%d", f->model.d[moNo]->s[sNo].yPosition) < 0) {
      GHMM_LOG(LERROR, "failed to write y position"); goto STOP;}
    if (xmlTextWriterEndElement(writer) < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement"); goto STOP;}
  }


  /* writing tied states */
  if (f->model.d[moNo]->model_type & GHMM_kTiedEmissions) {
    tied = f->model.d[moNo]->tied_to[sNo];
    if (tied != GHMM_kUntied) {
      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "tiedTo", "%d", tied);
      if (rc<0) {
	GHMM_LOG(LERROR, "failed to write tiedTo element");
	goto STOP;
      }
    }
  }


  return 0;
STOP:
  return -1;
#undef CUR_PROC
}

/* ========================================================================= */
static int writeContinuousStateContents(xmlTextWriterPtr writer, fileData_s * f,
				      int moNo, int sNo) {
#define CUR_PROC "writeContinuousStateContents"

  int i, rc;
  char * tmp;
  int allFixed = 0;

  /* writing continuous distribution */
  rc = xmlTextWriterStartElement(writer, BAD_CAST "mixture");
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
    goto STOP;
  }

  if (f->model.c[moNo]->s[sNo].fix)
    allFixed = 1;

  for(i=0;i<f->model.c[moNo]->s[sNo].M;i++){
    switch (f->model.c[moNo]->s[sNo].density[i]) {
      case normal:     
        rc = xmlTextWriterStartElement(writer, BAD_CAST "normal");
        if (rc < 0) {
          GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
          goto STOP;
        }
        if (writeDoubleAttribute(writer, "mean",f->model.c[moNo]->s[sNo].mue[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }
        if (writeDoubleAttribute(writer, "variance",f->model.c[moNo]->s[sNo].u[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }
        break;       
      case normal_left:     
        rc = xmlTextWriterStartElement(writer, BAD_CAST "normalTruncatedLeft");
        if (rc < 0) {
          GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
          goto STOP;
        }
        if (writeDoubleAttribute(writer, "mean",f->model.c[moNo]->s[sNo].mue[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }
        if (writeDoubleAttribute(writer, "variance",f->model.c[moNo]->s[sNo].u[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }
        if (writeDoubleAttribute(writer, "min",f->model.c[moNo]->s[sNo].a[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }    
        break;
      case normal_right:     
        rc = xmlTextWriterStartElement(writer, BAD_CAST "normalTruncatedRight");
        if (rc < 0) {
          GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
          goto STOP;
        }
        if (writeDoubleAttribute(writer, "mean",f->model.c[moNo]->s[sNo].mue[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }
        if (writeDoubleAttribute(writer, "variance",f->model.c[moNo]->s[sNo].u[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }
        if (writeDoubleAttribute(writer, "max",f->model.c[moNo]->s[sNo].a[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }    
        break;
      case uniform:     
        rc = xmlTextWriterStartElement(writer, BAD_CAST "uniform");
        if (rc < 0) {
          GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
          goto STOP;
        }
        if (writeDoubleAttribute(writer, "max",f->model.c[moNo]->s[sNo].mue[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }
        if (writeDoubleAttribute(writer, "min",f->model.c[moNo]->s[sNo].u[i])) {
  	  GHMM_LOG_QUEUED(LERROR); 
	  goto STOP;
        }    
        break;
      default:
        GHMM_LOG(LCRITIC, "invalid modelType");
        break;
    }
  
    /*optional values */ 
    if (allFixed || f->model.c[moNo]->s[sNo].mixture_fix[i]){
      if (writeIntAttribute(writer, "fixed", 1)) {
	GHMM_LOG_QUEUED(LERROR); 
	goto STOP;
      }
    }
    if (f->model.c[moNo]->s[sNo].M > 1){
      if (writeDoubleAttribute(writer, "prior",f->model.c[moNo]->s[sNo].c[i])) {
	GHMM_LOG_QUEUED(LERROR); 
	goto STOP;
      }
    }
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement");
      goto STOP;
    }

  }


  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement");
    goto STOP;
  }
  
  return 0;
STOP:
  return -1;
#undef CUR_PROC
}


/* ========================================================================= */
static int writeState(xmlTextWriterPtr writer, fileData_s * f, int moNo, int sNo) {
#define CUR_PROC "writeState"

  int rc;
  double w_pi;
  char * tmp, * w_desc=NULL;

  /* start state */
  rc = xmlTextWriterStartElement(writer, BAD_CAST "state");
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
    goto STOP;
  }

  /* write id attribute */
  if (writeIntAttribute(writer, "id", sNo))
    GHMM_LOG_QUEUED(LERROR);

  /* read state attribute from different model types */
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
  case GHMM_kDiscreteHMM:
    w_pi = f->model.d[moNo]->s[sNo].pi;
    /*w_desc = f->model.d[moNo]->s[sNo].desc;*/
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    /*
    w_pi = f->model.d[moNo]->s[sNo].pi;
    w_desc = f->model.d[moNo]->s[sNo];
    */
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    /*
    w_pi = f->model.d[moNo]->s[sNo].pi;
    w_desc = f->model.d[moNo]->s[sNo];
    */
    break;
  case GHMM_kContinuousHMM:
  case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):
    w_pi = f->model.c[moNo]->s[sNo].pi;
    /* w_desc = f->model.c[moNo]->s[sNo].desc; */
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");}

  /* write initial probability as attribute */
  if (writeDoubleAttribute(writer, "initial", w_pi))
    GHMM_LOG_QUEUED(LERROR);
  
  /* write state description */
  if (w_desc) {
    if (xmlTextWriterWriteAttribute(writer, BAD_CAST "desc", BAD_CAST w_desc))
      GHMM_LOG(LERROR, "writing state description failed");
  }


  /* write state contents for different model types */
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
  case GHMM_kDiscreteHMM:
    rc = writeDiscreteStateContents(writer, f, moNo, sNo);
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    /*
    rc = writeDiscreteStateContents(writer, f, moNo, sNo);
    */
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    /*
    rc = writeDiscreteStateContents(writer, f, moNo, sNo);
    */
    break;
  case GHMM_kContinuousHMM:
  case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):
    rc = writeContinuousStateContents(writer, f, moNo, sNo);
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");}

  if (rc) {
    GHMM_LOG(LERROR, "writing state contents failed");
    goto STOP;
  } 


  /* end state*/
  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement");
    goto STOP;
  }

  return 0;
STOP:
  return -1;
#undef CUR_PROC
}

/* ========================================================================= */
static int writeTransition(xmlTextWriterPtr writer, fileData_s * f, int moNo,
			   int sNo) {
#define CUR_PROC "writeTransition"

  int cos=1, i, j, rc; 
  int out_states, * out_id;
  double * * out_a;
  double * w_out_a;
  char * tmp;

  ARRAY_MALLOC(w_out_a, cos);

  /* write state contents for different model types */
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
  case GHMM_kDiscreteHMM:
    out_states = f->model.d[moNo]->s[sNo].out_states;
    out_id    = f->model.d[moNo]->s[sNo].out_id;
    out_a      = &(f->model.d[moNo]->s[sNo].out_a);
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    /*
    out_states = f->model.dp[moNo]->s[sNo].out_states;
    out_id    = f->model.dp[moNo]->s[sNo].out_id;
    out_a      = &(f->model.dp[moNo]->s[sNo].out_a);
    */
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    /*
    out_states = f->model.dp[moNo]->s[sNo].out_states;
    out_id    = f->model.dp[moNo]->s[sNo].out_id;
    out_a      = &(f->model.dp[moNo]->s[sNo].out_a);
    */
    break;
  case GHMM_kContinuousHMM:
  case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):
    out_states = f->model.c[moNo]->s[sNo].out_states;
    out_id    = f->model.c[moNo]->s[sNo].out_id;
    out_a      = f->model.c[moNo]->s[sNo].out_a;
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");}

  for (i=0; i<out_states; i++) {
    /* start state */
    rc = xmlTextWriterStartElement(writer, BAD_CAST "transition");
    if (rc < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
      goto STOP;
    }

    /* write source id (current state attribute */
    if (writeIntAttribute(writer, "source", sNo))
      GHMM_LOG_QUEUED(LERROR);

    /* write target id as attribute */
    if (writeIntAttribute(writer, "target", out_id[i]))
      GHMM_LOG_QUEUED(LERROR);

    for (j=0; j<cos; j++)
      w_out_a[j] = out_a[j][i];

    tmp = doubleArrayToCSV(w_out_a, cos);
    if (tmp) {
      rc = xmlTextWriterWriteElement(writer, BAD_CAST "probability", BAD_CAST tmp);
      m_free(tmp);
      if (rc<0) {
	GHMM_LOG(LERROR, "Error at xmlTextWriterWriteElement");
	goto STOP;
      }
    } else {
      GHMM_LOG(LERROR, "converting array to CSV failed");
      m_free(tmp);
      goto STOP;
    }

    /* end transition */
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement");
      goto STOP;
    }
  }

  return 0;
STOP:
  return -1;
#undef CUR_PROC
}


/* ========================================================================= */
static int writeHMM(xmlTextWriterPtr writer, fileData_s * f, int number) {
#define CUR_PROC "writeHMM"
  int rc, i, N;
  int w_cos;
  double w_prior;
  char * estr, * w_name, * w_type;

  rc = xmlTextWriterStartElement(writer, BAD_CAST "HMM");
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
    goto STOP;;
  }

  /* write HMM attributes applicable */
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
  case GHMM_kDiscreteHMM:
    w_name  = f->model.d[number]->name;
    w_type  = strModeltype(f->model.d[number]->model_type);
    w_prior = f->model.d[number]->prior;
    w_cos   = 1;
    N       =  f->model.d[number]->N;
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    /*
    w_name  = f->model.ds[number]->name;
    w_type  = strModeltype(f->model.ds[number]->model_type);
    w_prior = f->model.ds[number]->prior;
    w_cos   = 0;
    N       =  f->model.ds[number]->N;
    */
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    /*
    w_name  = f->model.dp[number]->name;
    w_type  = strModeltype(f->model.dp[number]->model_type);
    w_prior = f->model.dp[number]->prior;
    w_cos   = 0;
    N       =  f->model.dp[number]->N;
    */
    break;
  case GHMM_kContinuousHMM:
  case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):
    w_name  = NULL;/*f->model.c[number]->name;*/
    w_type  = strModeltype(f->modelType);
    w_prior = f->model.c[number]->prior;
    w_cos   = f->model.c[number]->cos;
    N       =  f->model.c[number]->N;
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");}

  if (w_name) {
    if (xmlTextWriterWriteAttribute(writer, BAD_CAST "name", BAD_CAST w_name))
      GHMM_LOG(LERROR, "writing HMM name failed");
  }
  if (xmlTextWriterWriteAttribute(writer, BAD_CAST "type", BAD_CAST w_type))
    GHMM_LOG(LERROR, "writing HMM type failed");

  if (w_prior >= 0.0)
    if (writeDoubleAttribute(writer, "prior", w_prior))
      GHMM_LOG_QUEUED(LERROR);
  
  if (w_cos > 1)
    if (writeIntAttribute(writer, "transitionClasses", w_cos))
      GHMM_LOG_QUEUED(LERROR);
  

  /* write alphabet if applicable */
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
  case GHMM_kDiscreteHMM:
    rc = writeAlphabet(writer, f->model.d[number]->alphabet);
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    /*rc = writeAlphabet(writer, f->model.ds[number]->alphabet);*/
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    /*rc = writeAlphabet(writer, f->model.dp[number]->alphabets[0]);
    if (rc) {
      GHMM_LOG(LERROR, "writing first alphabet of discrete pair HMM failed");
      goto STOP;
    }
    rc = writeAlphabet(writer, f->model.dp[number]->alphabets[1]);*/
    break;
  case GHMM_kContinuousHMM:
  case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):
    rc=0;
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");} 
  if (rc) {
    estr = ighmm_mprintf(NULL, 0, "writing alphabet for HMM %d (type %d) failed");
    GHMM_LOG(LERROR, estr);
    m_free(estr);
  }

  /* write background distributions if applicable */
  if ((f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
		       + GHMM_kPairHMM + GHMM_kContinuousHMM)) == GHMM_kDiscreteHMM) {
    if (writeBackground(writer, f->model.d[number]->bp))
      GHMM_LOG(LERROR, "writing of background distributions failed");
  }

  /* write all states */
  for (i=0; i<N; i++)
    if (writeState(writer, f, number, i)) {
      estr = ighmm_mprintf(NULL, 0, "writing of state %d in HMM %d failed", i, number);
      GHMM_LOG(LERROR, estr);
      m_free(estr);
      goto STOP;
    }

  /* write all outgoing transitions */
  for (i=0; i<N; i++)
    if (writeTransition(writer, f, number, i)) {
      estr = ighmm_mprintf(NULL, 0, "writing of state %d in HMM %d failed", i, number);
      GHMM_LOG(LERROR, estr);
      m_free(estr);
      goto STOP;
    }
      
  /*end HMM*/
  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement");
    goto STOP;
  }

  return 0;
STOP:
  return -1;
#undef CUR_PROC
}

/* ========================================================================= */
void writeHMMDocument(fileData_s * f, const char *file) {
#define CUR_PROC "writeHMMDocument"
  int rc, i;
  xmlTextWriterPtr writer;
  xmlChar *tmp;
  xmlDocPtr doc;
  
  /*
   * this initialize the library and check potential ABI mismatches
   * between the version it was compiled for and the actual shared
   * library used.
   */
  LIBXML_TEST_VERSION

  /* Create a new XmlWriter for DOM, with no compression. */
  writer = xmlNewTextWriterDoc(&doc, 0);
  if (writer == NULL) {
    GHMM_LOG(LERROR, "can not create the xml writer");
    goto STOP;
  }

  /* Start the document with the xml default for the version,
   * encoding ISO 8859-1 and the default for the standalone
   * declaration. */
  rc = xmlTextWriterStartDocument(writer, NULL, MY_ENCODING, NULL);
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartDocument\n");
    goto STOP;
  }

  rc = xmlTextWriterStartElement(writer, BAD_CAST "mixture");
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement\n");
    goto STOP;;
  }

  if (writeIntAttribute(writer, "noComponents", f->noModels)) {
    GHMM_LOG_QUEUED(LERROR);
    goto STOP;
  }

  for (i=0; i<f->noModels; i++)
    writeHMM(writer, f, i);

  /* Here we could close the elements ORDER and EXAMPLE using the
   * function xmlTextWriterEndElement, but since we do not want to
   * write any other elements, we simply call xmlTextWriterEndDocument,
   * which will do all the work. */
  rc = xmlTextWriterEndDocument(writer);
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndDocument");
    goto STOP;
  }

  xmlFreeTextWriter(writer);

  xmlSaveFormatFileEnc(file, doc, MY_ENCODING, 1);

STOP:
  xmlFreeDoc(doc);

  /*
   * Cleanup function for the XML library.
   */
  xmlCleanupParser();
  /*
   * this is to debug memory for regression tests
   */
  xmlMemoryDump();
#undef CUR_PROC
}

/**
 * ConvertInput:
 * @in: string in a given encoding
 * @encoding: the encoding used
 *
 * Converts @in into UTF-8 for processing with libxml2 APIs
 *
 * Returns the converted UTF-8 string, or NULL in case of error.
 */
static xmlChar *
ConvertInput(const char *in, const char *encoding)
{
    xmlChar *out;
    int ret;
    int size;
    int out_size;
    int temp;
    xmlCharEncodingHandlerPtr handler;

    if (in == 0)
        return 0;

    handler = xmlFindCharEncodingHandler(encoding);

    if (!handler) {
        printf("ConvertInput: no encoding handler found for '%s'\n",
               encoding ? encoding : "");
        return 0;
    }

    size = (int) strlen(in) + 1;
    out_size = size * 2 - 1;
    out = (unsigned char *) xmlMalloc((size_t) out_size);

    if (out != 0) {
        temp = size - 1;
        ret = handler->input(out, &out_size, (const xmlChar *) in, &temp);
        if ((ret < 0) || (temp - size + 1)) {
            if (ret < 0) {
                printf("ConvertInput: conversion wasn't successful.\n");
            } else {
                printf
                    ("ConvertInput: conversion wasn't successful. converted: %i octets.\n",
                     temp);
            }

            xmlFree(out);
            out = 0;
        } else {
            out = (unsigned char *) xmlRealloc(out, out_size + 1);
            out[out_size] = 0;  /*null terminating out */
        }
    } else {
        printf("ConvertInput: no mem\n");
    }

    return out;
}

#endif
