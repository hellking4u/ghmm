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

/* Bitmask to test the modeltype against to choose the type of the model pointer
   we use in the union */
#define PTR_TYPE_MASK (GHMM_kDiscreteHMM + GHMM_kTransitionClasses + GHMM_kPairHMM + GHMM_kContinuousHMM)


#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

#define MY_ENCODING "ISO-8859-1"

#define DTD_VERSION "1.0"

#define WRITE_DOUBLE_ATTRIBUTE(XMLW, NAME, VALUE)                             \
          if (0 > xmlTextWriterWriteFormatAttribute(XMLW, BAD_CAST (NAME),    \
                                                    "%.8f", (VALUE))) {       \
            estr=ighmm_mprintf(NULL,0, "failed to write attribute %s (%.8f)", \
                               (NAME), (VALUE));                              \
            GHMM_LOG(LERROR, estr); m_free(estr); goto STOP;} else



/* ========================================================================= */
static char * strModeltype(int modelType) {
#define CUR_PROC "strModelType"

  int end;
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

  /* overwrite the last space */
  end = strlen(mt);
  mt[end-1] = '\0';
  
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
    pos += sprintf(csv+pos, "%.8f, ", array[i]);
  }
  if (i<size-1) {
    GHMM_LOG(LERROR, "writing CSV failed");
    goto STOP;
  } else {
    pos += sprintf(csv+pos, "%.8f", array[i]);
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

  int i;
  char * estr;

  if (0 > xmlTextWriterStartElement(writer, BAD_CAST "alphabet")) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement");
    goto STOP;;
  }
  
  if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "id", "%d", alfa->id)) {
    estr = ighmm_mprintf(NULL, 0, "failed to write id-attribute for alphabet"
			 "with id %d", alfa->id);
    GHMM_LOG(LERROR, estr);
    m_free(estr);
  }
  
  for (i=0; i<alfa->size; i++) {
    if (0 > xmlTextWriterStartElement(writer, BAD_CAST "symbol")) {
      estr = ighmm_mprintf(NULL, 0, "failed to start symbol-tag no %d", i);
      GHMM_LOG(LERROR, estr);
      m_free(estr);
      goto STOP;
    }
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "code", "%d", i)) {
      estr = ighmm_mprintf(NULL, 0, "failed to write code-attribute for symbol %s"
			   "with code %d", alfa->symbols[i], i);
      GHMM_LOG(LERROR, estr);
      m_free(estr);
      goto STOP;
    }
    
    if (0 > xmlTextWriterWriteRaw(writer, BAD_CAST alfa->symbols[i])) {
      estr = ighmm_mprintf(NULL, 0, "failed to write symbol %s with code %d",
			   alfa->symbols[i], i);
      GHMM_LOG(LERROR, estr);
      m_free(estr);
      goto STOP;
    }
    
    if (0 > xmlTextWriterEndElement(writer)) {
      estr = ighmm_mprintf(NULL, 0, "failed to end symbol-tag no %d", i);
      GHMM_LOG(LERROR, estr);
      m_free(estr);
      goto STOP;
    }
  }

  if (0 > xmlTextWriterEndElement(writer)) {
    GHMM_LOG(LERROR, "Error at ending alphabet");
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
  
  int i;
  char * tmp=NULL;
  char * estr;

  for (i=0; i<bg->n; i++) {

    if (0 > xmlTextWriterStartElement(writer, BAD_CAST "background")) {
      estr = ighmm_mprintf(NULL, 0, "Error at starting backgroung %d", i);
      GHMM_LOG(LERROR, estr);
      m_free(estr);
      return -1;
    }
    
    if (0 > xmlTextWriterWriteAttribute(writer, BAD_CAST "key", BAD_CAST bg->name[i]))
      GHMM_LOG(LERROR, "Error at writing background key");

    if (0 < bg->order[i])
      if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "order", "%d", bg->order[i]))
	GHMM_LOG(LERROR, "can't write background order attribute");
    
    tmp = doubleArrayToCSV(bg->b[i], pow(bg->m, bg->order[i]+1));
    if (tmp) {
      if (0 > xmlTextWriterWriteRaw(writer, BAD_CAST tmp)) {
	GHMM_LOG(LERROR, "Error at xmlTextWriterWriteRaw while writing"
		 "background distribution CSV");
	m_free(tmp);
	return -1;
      }
      m_free(tmp);
    } else {
      GHMM_LOG(LERROR, "converting array to CSV failed for background distribution");
      return -1;
    }

    if (0 > xmlTextWriterEndElement(writer)) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement while ending"
	       "background distribution");
      return -1;
    }
  }
  return 0;
#undef CUR_PROC
}

/* ========================================================================= */
static int writeDiscreteStateContents(xmlTextWriterPtr writer, fileData_s * f,
				      int moNo, int sNo) {
#define CUR_PROC "writeDiscreteStateContents"

  int bgId, cLabel, rc, order, tied;
  char * tmp=NULL;

  /* writing discrete distribution */
  if (0 > xmlTextWriterStartElement(writer, BAD_CAST "discrete")) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (discrete)");
    goto STOP;
  }

  if (0 > xmlTextWriterWriteAttribute(writer, BAD_CAST "id", BAD_CAST "0")) {
    GHMM_LOG(LERROR, "failed to write alphabet id");
    goto STOP;
  }

  if (f->model.d[moNo]->s[sNo].fix)
    if (0 > xmlTextWriterWriteAttribute(writer, BAD_CAST "fixed", BAD_CAST "1")) {
      GHMM_LOG(LERROR, "failed to write fixed attriute");
      goto STOP;
    }

  if ((f->model.d[moNo]->model_type & GHMM_kHigherOrderEmissions)
      && f->model.d[moNo]->order[sNo]) {
    order = f->model.d[moNo]->order[sNo];
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "order", "%d", order)) {
      GHMM_LOG(LERROR, "failed to write order attribute for discrete distribution");
      goto STOP;
    }
  } else
    order = 0;

  tmp = doubleArrayToCSV(f->model.d[moNo]->s[sNo].b, pow(f->model.d[moNo]->M, order+1));
  if (tmp) {
    if (0 > xmlTextWriterWriteRaw(writer, BAD_CAST tmp)) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterWriteRaw while writing"
	       "discrete distribution CSV");
      m_free(tmp);
      goto STOP;
    }
    m_free(tmp);
  } else {
    GHMM_LOG(LERROR, "converting array to CSV failed for discrete distribution");
    goto STOP;
  }
  
  /* end discrete distribution */
  if (0 > xmlTextWriterEndElement(writer)) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (discrete)");
    goto STOP;
  }

  /* writing backgroung key */
  if (f->model.d[moNo]->model_type & GHMM_kBackgroundDistributions) {
    bgId = f->model.d[moNo]->background_id[sNo];
    if (bgId != GHMM_kNoBackgroundDistribution) {
      if (f->model.d[moNo]->bp->name[bgId]) {
	rc = xmlTextWriterWriteElement(writer, BAD_CAST "backgroundKey",
				       BAD_CAST f->model.d[moNo]->bp->name[bgId]);
	if (rc<0) {
	  GHMM_LOG(LERROR, "Error at xmlTextWriterWriteElement (backgroundKey)");
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
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "class", "%d", cLabel);
    if (rc<0) {
      GHMM_LOG(LERROR, "failed to write class label");
      goto STOP;
    }
  }
  
  /* duration (not implemented yet, maybe never */
#if 0
  if (f->model.d[moNo]->model_type & GHMM_kDurations) {
    if (f->model.d[moNo]->duration[sNo] > 0) {
      rc = xmlTextWriterWriteElement(writer, BAD_CAST "duration",
				     BAD_CAST f->model.d[moNo]->duration[sNo]);
      if (rc<0) {
	GHMM_LOG(LERROR, "Error at xmlTextWriterWriteElement (duration)");
	goto STOP;
      }
    }
  }
#endif

  /* writing positions */
  if ((f->model.d[moNo]->s[sNo].xPosition > 0)
      && (f->model.d[moNo]->s[sNo].xPosition > 0)) {
    if (xmlTextWriterStartElement(writer, BAD_CAST "position") < 0) {
      GHMM_LOG(LERROR, "failed to start position element (position)"); goto STOP;}
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "x", "%d",
					  f->model.d[moNo]->s[sNo].xPosition)) {
      GHMM_LOG(LERROR, "failed to write x position"); goto STOP;}
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "y", "%d",
					  f->model.d[moNo]->s[sNo].yPosition)) {
      GHMM_LOG(LERROR, "failed to write y position"); goto STOP;}
    if (xmlTextWriterEndElement(writer) < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (position)"); goto STOP;}
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
static int writeDiscreteSwitchingStateContents(xmlTextWriterPtr writer,
					       fileData_s * f, int moNo,
					       int sNo) {
#define CUR_PROC "writeDiscreteSwitchingStateContents"

  int bgId, cLabel, rc, order, tied;
  char * tmp=NULL;

  /* writing discrete distribution */
  if (0 > xmlTextWriterStartElement(writer, BAD_CAST "discrete")) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (discrete)");
    goto STOP;
  }

  if (0 > xmlTextWriterWriteAttribute(writer, BAD_CAST "id", BAD_CAST "0")) {
    GHMM_LOG(LERROR, "failed to write alphabet id");
    goto STOP;
  }

  if (f->model.ds[moNo]->s[sNo].fix)
    if (0 > xmlTextWriterWriteAttribute(writer, BAD_CAST "fixed", BAD_CAST "1")) {
      GHMM_LOG(LERROR, "failed to write fixed attriute");
      goto STOP;
    }

  if ((f->model.ds[moNo]->model_type & GHMM_kHigherOrderEmissions)
      && f->model.ds[moNo]->order[sNo]) {
    order = f->model.ds[moNo]->order[sNo];
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "order", "%d", order)) {
      GHMM_LOG(LERROR, "failed to write order attribute for discrete distribution");
      goto STOP;
    }
  } else
    order = 0;

  tmp = doubleArrayToCSV(f->model.ds[moNo]->s[sNo].b, pow(f->model.ds[moNo]->M, order+1));
  if (tmp) {
    if (0 > xmlTextWriterWriteRaw(writer, BAD_CAST tmp)) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterWriteRaw while writing"
	       "discrete distribution CSV");
      m_free(tmp);
      goto STOP;
    }
    m_free(tmp);
  } else {
    GHMM_LOG(LERROR, "converting array to CSV failed for discrete distribution");
    goto STOP;
  }
  
  /* end discrete distribution */
  if (0 > xmlTextWriterEndElement(writer)) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (discrete)");
    goto STOP;
  }

  /* writing backgroung key */
  if (f->model.ds[moNo]->model_type & GHMM_kBackgroundDistributions) {
    bgId = f->model.ds[moNo]->background_id[sNo];
    if (bgId != GHMM_kNoBackgroundDistribution) {
      if (f->model.ds[moNo]->bp->name[bgId]) {
	rc = xmlTextWriterWriteElement(writer, BAD_CAST "backgroundKey",
				       BAD_CAST f->model.ds[moNo]->bp->name[bgId]);
	if (rc<0) {
	  GHMM_LOG(LERROR, "Error at xmlTextWriterWriteElement (backgroundKey)");
	  goto STOP;
	}
      } else {
	GHMM_LOG(LERROR, "background name is NULL pointer, invalid model");
	goto STOP;
      }
    }
  }

  /* writing class label */
  if (f->model.ds[moNo]->model_type & GHMM_kLabeledStates) {
    cLabel = f->model.ds[moNo]->label[sNo];
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "class", "%d", cLabel);
    if (rc<0) {
      GHMM_LOG(LERROR, "failed to write class label");
      goto STOP;
    }
  }
  
  /* duration (not implemented yet, maybe never */
#if 0
  if (f->model.ds[moNo]->model_type & GHMM_kDurations) {
    if (f->model.ds[moNo]->duration[sNo] > 0) {
      rc = xmlTextWriterWriteElement(writer, BAD_CAST "duration",
				     BAD_CAST f->model.ds[moNo]->duration[sNo]);
      if (rc<0) {
	GHMM_LOG(LERROR, "Error at xmlTextWriterWriteElement (duration)");
	goto STOP;
      }
    }
  }
#endif

  /* writing positions */
  if ((f->model.ds[moNo]->s[sNo].xPosition > 0)
      && (f->model.ds[moNo]->s[sNo].xPosition > 0)) {
    if (xmlTextWriterStartElement(writer, BAD_CAST "position") < 0) {
      GHMM_LOG(LERROR, "failed to start position element (position)"); goto STOP;}
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "x", "%d",
					  f->model.ds[moNo]->s[sNo].xPosition)) {
      GHMM_LOG(LERROR, "failed to write x position"); goto STOP;}
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "y", "%d",
					  f->model.ds[moNo]->s[sNo].yPosition)) {
      GHMM_LOG(LERROR, "failed to write y position"); goto STOP;}
    if (xmlTextWriterEndElement(writer) < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (position)"); goto STOP;}
  }


  /* writing tied states */
  if (f->model.ds[moNo]->model_type & GHMM_kTiedEmissions) {
    tied = f->model.ds[moNo]->tied_to[sNo];
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

  int i;
  int allFixed = 0;
  char * estr; /* needed for WRITE_DOUBLE_ATTRIBUTE macro */

  /* writing continuous distribution */
  if (0 > xmlTextWriterStartElement(writer, BAD_CAST "mixture")) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (mixture)");
    goto STOP;
  }

  if (f->model.c[moNo]->s[sNo].fix)
    allFixed = 1;

  for(i=0;i<f->model.c[moNo]->s[sNo].M;i++){
    switch (f->model.c[moNo]->s[sNo].density[i]) {
      case normal:
        if (0 > xmlTextWriterStartElement(writer, BAD_CAST "normal")) {
          GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (normal)");
          goto STOP;
        }
	WRITE_DOUBLE_ATTRIBUTE(writer, "mean", f->model.c[moNo]->s[sNo].mue[i]);
	WRITE_DOUBLE_ATTRIBUTE(writer, "variance", f->model.c[moNo]->s[sNo].u[i]);
        break;       
      case normal_left:     
        if (0 > xmlTextWriterStartElement(writer, BAD_CAST "normalTruncatedLeft")) {
          GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (normalTruncatedLeft)");
          goto STOP;
        }
	WRITE_DOUBLE_ATTRIBUTE(writer, "mean", f->model.c[moNo]->s[sNo].mue[i]);
	WRITE_DOUBLE_ATTRIBUTE(writer, "variance", f->model.c[moNo]->s[sNo].u[i]);
        WRITE_DOUBLE_ATTRIBUTE(writer, "min", f->model.c[moNo]->s[sNo].a[i]);
        break;
      case normal_right:     
        if (0 > xmlTextWriterStartElement(writer, BAD_CAST "normalTruncatedRight")) {
          GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (normalTruncatedRight)");
          goto STOP;
        }
	WRITE_DOUBLE_ATTRIBUTE(writer, "mean", f->model.c[moNo]->s[sNo].mue[i]);
	WRITE_DOUBLE_ATTRIBUTE(writer, "variance", f->model.c[moNo]->s[sNo].u[i]);
        WRITE_DOUBLE_ATTRIBUTE(writer, "max", f->model.c[moNo]->s[sNo].a[i]);
        break;
      case uniform:
        if (0 > xmlTextWriterStartElement(writer, BAD_CAST "uniform")) {
          GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (uniform)");
          goto STOP;
        }
        WRITE_DOUBLE_ATTRIBUTE(writer, "min", f->model.c[moNo]->s[sNo].u[i]);
        WRITE_DOUBLE_ATTRIBUTE(writer, "max", f->model.c[moNo]->s[sNo].mue[i]);
        break;
      default:
        GHMM_LOG(LERROR, "invalid density");
	goto STOP;
    }
  
    /*optional values */ 
    if (allFixed || f->model.c[moNo]->s[sNo].mixture_fix[i]){
      if (0 > xmlTextWriterWriteAttribute(writer, BAD_CAST "fixed", BAD_CAST "1")) {
	GHMM_LOG(LERROR, "failed to set fixed attribute"); 
	goto STOP;
      }
    }
    if (f->model.c[moNo]->s[sNo].M > 1){
      WRITE_DOUBLE_ATTRIBUTE(writer, "prior", f->model.c[moNo]->s[sNo].c[i]);
    }
    if (0 > xmlTextWriterEndElement(writer)) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (all densities)");
      goto STOP;
    }
  }

  /* end mixture tag */
  if (0 > xmlTextWriterEndElement(writer)) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (mixture)");
    goto STOP;
  }
  
  /* writing positions */
  if ((f->model.c[moNo]->s[sNo].xPosition > 0)
      && (f->model.c[moNo]->s[sNo].yPosition > 0)) {
    if (xmlTextWriterStartElement(writer, BAD_CAST "position") < 0) {
      GHMM_LOG(LERROR, "failed to start position element (position)"); goto STOP;}
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "x", "%d",
					  f->model.c[moNo]->s[sNo].xPosition)) {
      GHMM_LOG(LERROR, "failed to write x position"); goto STOP;    }
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "y", "%d",
					  f->model.c[moNo]->s[sNo].yPosition)) {
      GHMM_LOG(LERROR, "failed to write y position"); goto STOP;}
    if (xmlTextWriterEndElement(writer) < 0) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (position)"); goto STOP;}
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
  unsigned char * w_desc=NULL;
  char * estr; /* needed for WRITE_DOUBLE_ATTRIBUTE macro */

  /* start state */
  if (0 > xmlTextWriterStartElement(writer, BAD_CAST "state")) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (state)");
    goto STOP;
  }

  /* write id attribute */
  if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "id", "%d", sNo))
    GHMM_LOG(LERROR, "failed to write statte id attribute");

  /* read state attribute from different model types */
  switch (f->modelType & PTR_TYPE_MASK) {
  case GHMM_kDiscreteHMM:
    w_pi = f->model.d[moNo]->s[sNo].pi;
    w_desc = f->model.d[moNo]->s[sNo].desc;
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    w_pi = f->model.ds[moNo]->s[sNo].pi;
    w_desc = f->model.ds[moNo]->s[sNo].desc;
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
    w_desc = f->model.c[moNo]->s[sNo].desc;
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");}

  /* write initial probability as attribute */
  WRITE_DOUBLE_ATTRIBUTE(writer, "initial", w_pi);
  
  /* write state description */
  if (w_desc) {
    if (xmlTextWriterWriteAttribute(writer, BAD_CAST "desc", w_desc))
      GHMM_LOG(LERROR, "writing state description failed");
  }

  /* write state contents for different model types */
  switch (f->modelType & PTR_TYPE_MASK) {
  case GHMM_kDiscreteHMM:
    rc = writeDiscreteStateContents(writer, f, moNo, sNo);
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    rc = writeDiscreteSwitchingStateContents(writer, f, moNo, sNo);
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    /*
    rc = writeDiscretePairStateContents(writer, f, moNo, sNo);
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
  if (0 > xmlTextWriterEndElement(writer)) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (state)");
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

  int cos, i, j, rc=0; 
  int out_states, * out_id;
  double * * out_a;
  double * w_out_a;
  char * tmp;

  ARRAY_MALLOC(w_out_a, cos);

  /* write state contents for different model types */
  switch (f->modelType & PTR_TYPE_MASK) {
  case GHMM_kDiscreteHMM:
    out_states = f->model.d[moNo]->s[sNo].out_states;
    out_id     = f->model.d[moNo]->s[sNo].out_id;
    out_a      = &(f->model.d[moNo]->s[sNo].out_a);
    cos        = 1;
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    out_states = f->model.ds[moNo]->s[sNo].out_states;
    out_id     = f->model.ds[moNo]->s[sNo].out_id;
    out_a      = f->model.ds[moNo]->s[sNo].out_a;
    cos        = f->model.ds[moNo]->cos;
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    /*
    out_states = f->model.dp[moNo]->s[sNo].out_states;
    out_id     = f->model.dp[moNo]->s[sNo].out_id;
    out_a      = f->model.dp[moNo]->s[sNo].out_a;
    cos        = f->model.dp[moNo]->cos;
    */
    break;
  case GHMM_kContinuousHMM:
  case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):
    out_states = f->model.c[moNo]->s[sNo].out_states;
    out_id     = f->model.c[moNo]->s[sNo].out_id;
    out_a      = f->model.c[moNo]->s[sNo].out_a;
    cos        = f->model.c[moNo]->cos;
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");}

  for (i=0; i<out_states; i++) {
    if (0 > xmlTextWriterStartElement(writer, BAD_CAST "transition")) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (transition)");
      goto STOP;
    }

    /* write source id (current state attribute */
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "source", "%d", sNo))
      GHMM_LOG(LERROR, "failed to write transition source attribute");

    /* write target id as attribute */
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "target", "%d", out_id[i]))
      GHMM_LOG(LERROR, "failed to write transition target attribute");

    for (j=0; j<cos; j++)
      w_out_a[j] = out_a[j][i];

    tmp = doubleArrayToCSV(w_out_a, cos);
    if (tmp) {
      if (0 > xmlTextWriterWriteElement(writer, BAD_CAST "probability", BAD_CAST tmp)) {
	GHMM_LOG(LERROR, "Error at xmlTextWriterWriteElement (transition probabilities)");
	m_free(tmp);
	goto STOP;
      }
      m_free(tmp);      
    } else {
      GHMM_LOG(LERROR, "converting transition probabilities array to CSV failed");
      goto STOP;
    }

    /* end transition */
    if (rc < 0 > xmlTextWriterEndElement(writer)) {
      GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (transition)");
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
  int rc=0, i, N;
  int w_cos;
  double w_prior;
  unsigned char * w_name;
  char * w_type;
  char * estr; /* needed for WRITE_DOUBLE_ATTRIBUTE macro */

  /* start HMM */
  if (0 > xmlTextWriterStartElement(writer, BAD_CAST "HMM")) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (HMM)");
    goto STOP;;
  }

  /* write HMM attributes applicable */
  switch (f->modelType & PTR_TYPE_MASK) {
  case GHMM_kDiscreteHMM:
    w_name  = f->model.d[number]->name;
    w_type  = strModeltype(f->model.d[number]->model_type);
    w_prior = f->model.d[number]->prior;
    N       = f->model.d[number]->N;
    w_cos   = 1;
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    w_name  = f->model.ds[number]->name;
    w_type  = strModeltype(f->model.ds[number]->model_type);
    w_prior = f->model.ds[number]->prior;
    N       = f->model.ds[number]->N;
    w_cos   = 0;
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    /*
    w_name  = f->model.dp[number]->name;
    w_type  = strModeltype(f->model.dp[number]->model_type);
    w_prior = f->model.dp[number]->prior;
    N       = f->model.dp[number]->N;
    w_cos   = 0;
    */
    break;
  case GHMM_kContinuousHMM:
  case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):
    w_name  = f->model.c[number]->name;
    w_type  = strModeltype(f->modelType);
    w_prior = f->model.c[number]->prior;
    N       = f->model.c[number]->N;
    w_cos   = f->model.c[number]->cos;
    break;
  default:
    GHMM_LOG(LERROR, "invalid modelType");
    goto STOP;}

  if (w_name) {
    if (xmlTextWriterWriteAttribute(writer, BAD_CAST "name", w_name))
      GHMM_LOG(LERROR, "writing HMM name failed");
  }
  if (xmlTextWriterWriteAttribute(writer, BAD_CAST "type", BAD_CAST w_type))
    GHMM_LOG(LERROR, "writing HMM type failed");

  if (w_prior >= 0.0) {
    WRITE_DOUBLE_ATTRIBUTE(writer, "prior", w_prior);
  }

  if (w_cos > 1)
    if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "transitionClasses",
					      "%d", w_cos))
      GHMM_LOG(LERROR, "failed to write no of transitionClasses");
  

  /* write alphabet if applicable */
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM)) {
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
  }

  if (rc) {
    estr = ighmm_mprintf(NULL, 0, "writing alphabet for HMM %d (type %d) failed");
    GHMM_LOG(LERROR, estr);
    m_free(estr);
  }

  /* write background distributions if applicable */
  if ((f->modelType & PTR_TYPE_MASK) == GHMM_kDiscreteHMM
      && f->modelType & GHMM_kBackgroundDistributions) {
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
      estr = ighmm_mprintf(NULL, 0, "writing transitions of state %d in HMM %d failed",
			   i, number);
      GHMM_LOG(LERROR, estr);
      m_free(estr);
      goto STOP;
    }
      
  /*end HMM*/
  if (0 > xmlTextWriterEndElement(writer)) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndElement (HMM)");
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

  /* indenting writer to circumvent no space between SYSTEM and PUBLIC identifier */
  xmlTextWriterSetIndent(writer, 1);

  /* Start the document with the xml default for the version,
   * encoding ISO 8859-1 and the default for the standalone
   * declaration. */
  rc = xmlTextWriterStartDocument(writer, NULL, MY_ENCODING, NULL);
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartDocument\n");
    goto STOP;
  }

  /* Set the Document type declaration at the beginning of the document */
  rc = xmlTextWriterWriteDTD(writer, BAD_CAST "mixture",
			     BAD_CAST "-//ghmm.org//DOCUMENT ghmm V"DTD_VERSION"//EN",
			     BAD_CAST "http://ghmm.sourceforge.net/xml/"DTD_VERSION"/ghmm.dtd",
			     NULL);
  if (rc < 0) {
    GHMM_LOG(LERROR, "failed to write the DocType"); goto STOP;}

  /* start real contents */
  if (0 > xmlTextWriterStartElement(writer, BAD_CAST "mixture")) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterStartElement (mixture)");
    goto STOP;;
  }

  if (xmlTextWriterWriteAttribute(writer, BAD_CAST "version", BAD_CAST DTD_VERSION) < 0) {
    GHMM_LOG(LERROR, "failed to write version 1.0"); goto STOP;}

  if (0 > xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "noComponents", "%d", f->noModels)) {
    GHMM_LOG(LERROR, "failed to write the number of components"); goto STOP;}

  /* write all models */
  for (i=0; i<f->noModels; i++)
    writeHMM(writer, f, i);

  /* end mixture */
  if (0 > xmlTextWriterEndDocument(writer)) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndDocument (mixture)");
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

#endif
