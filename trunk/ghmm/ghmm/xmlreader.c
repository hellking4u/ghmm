/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/xmlreader.c
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

#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>

#include <libxml/xmlmemory.h>
#include <libxml/tree.h>
#include <libxml/parser.h>

#include "ghmm.h"
#include "ghmm_internals.h"
#include "mes.h"
#include "mprintf.h"
#include "xmlreader.h"

/* we should not need more than to alphabets, no plan to implement triple HMMs */
#define MAX_ALPHABETS 2 


/* holds all valid modeltypes sorted */
static int validModelTypes[28] = {
  (GHMM_kDiscreteHMM),
  (GHMM_kDiscreteHMM + GHMM_kLeftRight),
  (GHMM_kDiscreteHMM + GHMM_kSilentStates),
  (GHMM_kDiscreteHMM + GHMM_kTiedEmissions),
  (GHMM_kDiscreteHMM + GHMM_kTiedEmissions + GHMM_kSilentStates),
  (GHMM_kDiscreteHMM + GHMM_kHigherOrderEmissions),
  (GHMM_kDiscreteHMM + GHMM_kHigherOrderEmissions + GHMM_kSilentStates),
  (GHMM_kDiscreteHMM + GHMM_kHigherOrderEmissions + GHMM_kTiedEmissions),
  (GHMM_kDiscreteHMM + GHMM_kHigherOrderEmissions + GHMM_kTiedEmissions	+ GHMM_kSilentStates),
  (GHMM_kDiscreteHMM + GHMM_kBackgroundDistributions),
  (GHMM_kDiscreteHMM + GHMM_kBackgroundDistributions + GHMM_kSilentStates),
  (GHMM_kDiscreteHMM + GHMM_kBackgroundDistributions + GHMM_kTiedEmissions),
  (GHMM_kDiscreteHMM + GHMM_kBackgroundDistributions + GHMM_kTiedEmissions + GHMM_kSilentStates),
  (GHMM_kDiscreteHMM + GHMM_kBackgroundDistributions + GHMM_kHigherOrderEmissions + GHMM_kSilentStates),
  (GHMM_kDiscreteHMM + GHMM_kBackgroundDistributions + GHMM_kHigherOrderEmissions + GHMM_kTiedEmissions),
  (GHMM_kDiscreteHMM + GHMM_kBackgroundDistributions + GHMM_kHigherOrderEmissions + GHMM_kTiedEmissions + GHMM_kSilentStates),
  (GHMM_kDiscreteHMM + GHMM_kLabeledStates),
  (GHMM_kDiscreteHMM + GHMM_kLabeledStates + GHMM_kTiedEmissions),
  (GHMM_kDiscreteHMM + GHMM_kLabeledStates + GHMM_kHigherOrderEmissions),
  (GHMM_kDiscreteHMM + GHMM_kLabeledStates + GHMM_kHigherOrderEmissions + GHMM_kTiedEmissions),
  (GHMM_kDiscreteHMM + GHMM_kLabeledStates + GHMM_kBackgroundDistributions),
  (GHMM_kDiscreteHMM + GHMM_kLabeledStates + GHMM_kBackgroundDistributions + GHMM_kTiedEmissions),
  (GHMM_kDiscreteHMM + GHMM_kLabeledStates + GHMM_kBackgroundDistributions + GHMM_kHigherOrderEmissions + GHMM_kTiedEmissions),
  (GHMM_kDiscreteHMM + GHMM_kTransitionClasses),
  (GHMM_kContinuousHMM),
  (GHMM_kContinuousHMM + GHMM_kTransitionClasses),
  (GHMM_kPairHMM + GHMM_kDiscreteHMM),
  (GHMM_kPairHMM + GHMM_kDiscreteHMM + GHMM_kTransitionClasses)
};


/*===========================================================================*/
static int getIntAttribute(xmlNodePtr node, const xmlChar *name, int *error) {
  xmlChar *attr;
  int value = -3894;

  if ((attr = xmlGetProp(node, name)) != NULL) {
    value = atoi((char*)attr);
    xmlFree(attr);
    *error = 0;
  } else {
    *error = 1;
  }
  return value;
}

/*===========================================================================*/
static double getDoubleAttribute(xmlNodePtr node, const xmlChar *name,
				 int *error) {
  xmlChar *attr;
  double value = 0.0;

  if ((attr = xmlGetProp(node, name)) != NULL) {
    value = atof((char*)attr);
    xmlFree(attr);
    *error = 0;
  } else {
    *error = 1;
  }
  return value;
}

/*===========================================================================*/
/* Caller owns return value */
static xmlChar* getXMLCharAttribute(xmlNodePtr node, const xmlChar *name,
				    int *error) {

  xmlChar *attr;

  if ((attr = xmlGetProp(node, name)) != NULL) {
    *error = 0;
    return attr;
  } else {
    *error = 1;
    return NULL;
  }
}


/*===========================================================================*/
static int parseCSVList(const char * data, unsigned int size, double * array) {
#define CUR_PROC "parseCSVList"

  int retval=0;
  int i;
  char * * next, * estr;

  ARRAY_CALLOC(next, 1);

  for (i=0; i<size; i++) {
    array[i] = strtod(data, next);
    if (data == *next) {
      estr = ighmm_mprintf(NULL, 0, "error in parsing CSV. entry %d of %d. (%s)", i, size, *next);
      GHMM_LOG(LERROR, estr);
      m_free(estr);
      retval=-1;
      break;
    }
    if (next)
      data = *next+1;
    else
      break;
  }

  if (i != size) {
    retval=-1;
    estr = ighmm_mprintf(NULL, 0, "error in parsing CSV. sizes do not match (%d != %d)", i, size);
    GHMM_LOG(LERROR, estr);
    m_free(estr);
    
  }

STOP:
  m_free(next);
  return retval;
#undef CUR_PROC
}

/*===========================================================================*/
static int matchModelType(const char * data, unsigned int size) {
#define CUR_PROC "matchModelType"

  if (!strncmp(data, "left-right", size))
    return GHMM_kLeftRight;

  if (!strncmp(data, "silent", size))
    return GHMM_kSilentStates;

  if (!strncmp(data, "tied", size))
    return GHMM_kTiedEmissions;

  if (!strncmp(data, "higher-order", size))
    return GHMM_kHigherOrderEmissions;

  if (!strncmp(data, "background", size))
    return GHMM_kBackgroundDistributions;

  if (!strncmp(data, "labeled", size))
    return GHMM_kLabeledStates;

  if (!strncmp(data, "transition-classes", size))
    return GHMM_kTransitionClasses;

  if (!strncmp(data, "discrete", size))
    return GHMM_kDiscreteHMM;

  if (!strncmp(data, "continuous", size))
    return GHMM_kContinuousHMM;

  if (!strncmp(data, "pair", size))
    return GHMM_kPairHMM;

  return INT_MIN;


#undef CUR_PROC
}
/*===========================================================================*/
static int parseModelType(const char * data, unsigned int size) {
#define CUR_PROC "parseModelType"

  int i, modelType=0;
  const char * end = data;
  char * str;

  while ((end = strchr(data, ' '))) {
    modelType += matchModelType(data, end-data);
    size -= (end-data)+1;
    data = end+1;
  }
  
  modelType += matchModelType(data, size);

  for (i=0; i<sizeof(validModelTypes); i++) {
    if (modelType == validModelTypes[i])
      break;
  }
  if (i == sizeof(validModelTypes)) {
    str = ighmm_mprintf(NULL, 0, "%d is no known valid model type", modelType);
    GHMM_LOG(LERROR, str);
    m_free(str);
    return -1;
  }

  return modelType;
#undef CUR_PROC
}


/*===========================================================================*/
static alphabet_s * parseAlphabet(xmlDocPtr doc, xmlNodePtr cur, fileData_s * f) {
#define CUR_PROC "parseAlphabet"
  
  char * str;
  int M, code, error;

  xmlNodePtr symbol;
  xmlChar * s;
  alphabet_s * alfa;

  ARRAY_CALLOC(alfa, 1);

  symbol = cur->children;
  M=0;
  while (symbol!=NULL) {
    if ((!xmlStrcmp(symbol->name, (const xmlChar *)"symbol"))) {
      code = getIntAttribute(symbol, (const xmlChar *)"code", &error);
      if (error || code!=M) {
	str = ighmm_mprintf(NULL, 0, "non consecutive code %d == %d", code, M);
	GHMM_LOG(LERROR, str);
	m_free(str);
	goto STOP;
      } else
	M++;
    }
    symbol=symbol->next;
  }

  alfa->size = M;
  ARRAY_MALLOC(alfa->symbols, M);

  symbol = cur->xmlChildrenNode;
  M=0;
  while (symbol!=NULL) {
    if ((!xmlStrcmp(symbol->name, (const xmlChar *)"symbol"))) {
      s = xmlNodeGetContent(cur);
      alfa->symbols[M++] = (char *)s;
    }
    symbol=symbol->next;
  }

  return alfa;
STOP:
  m_free(alfa->symbols);
  m_free(alfa)
  return NULL;
#undef CUR_PROC
}

/*===========================================================================*/
static int parseBackground(xmlDocPtr doc, xmlNodePtr cur, fileData_s * f, int modelNo) {
#define CUR_PROC "parseBackground"

  int error, order;
  int bgNr;
  double * b;
  char * s;

  assert(f->modelType & GHMM_kDiscreteHMM);

  bgNr = f->model.d[modelNo]->bp->n++;

  /* get order */
  order = getIntAttribute(cur, (const xmlChar *)"order", &error);
  if (error)
    order=0;
  else if (order && !(f->modelType & GHMM_kHigherOrderEmissions)) {
    GHMM_LOG(LERROR, "background distribution has order > 0, but model is not higher order");
    goto STOP;
  }
  f->model.d[modelNo]->bp->order[bgNr] = order;

  /* get name */
  s = (char *)getXMLCharAttribute(cur, (const xmlChar *)"name", &error);
  f->model.d[modelNo]->bp->name[bgNr] = s;

  /* get distribution */
  s = (char *)xmlNodeGetContent(cur);

  ARRAY_MALLOC(b, pow(f->model.d[modelNo]->bp->m, order+1));
  if (-1 !=  parseCSVList(s, pow(f->model.d[modelNo]->bp->m, order+1), b))
    f->model.d[modelNo]->bp->b[bgNr] = b;
  else {
    GHMM_LOG(LERROR, "Can not parse background CSV list.");
    goto STOP;
  }

  return 0;
STOP:
  m_free(b);
  return -1;
#undef CUR_PROC
}

/*===========================================================================*/
static int parseState(xmlDocPtr doc, xmlNodePtr cur, fileData_s * f, int * inDegree, int * outDegree, int modelNo) {
#define CUR_PROC "parseState"

  int i, error, order=-28, state=-1442, fixed=-985, tied=-9354, M, aprox;
  double pi, prior;
  double * emissions;
  char * desc, * s, * estr, * * serror;


  xmlNodePtr elem, child;

  state = getIntAttribute(cur, (const xmlChar *)"id", &error);
  pi = getDoubleAttribute(cur, (const xmlChar *)"initial", &error);
  if (error) {
    estr = ighmm_mprintf(NULL, 0, "can't read required intial probability for"
			 "state %d", state);
    GHMM_LOG(LERROR, estr);
    goto STOP;
  } else
    
  desc = (char *)getXMLCharAttribute(cur, (const xmlChar *)"desc", &error);

  elem = cur->children;
  while (elem!=NULL) {
    /* ======== silent state ============================================== */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"silent"))) {
      if (f->modelType & GHMM_kDiscreteHMM)
	f->model.d[modelNo]->silent[state] = 1;
    }

    /* ======== discrete state (possible higher order) ==================== */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"discrete"))) {
      assert((f->modelType & GHMM_kDiscreteHMM) && ((f->modelType & GHMM_kPairHMM) == 0));
      
      f->model.d[modelNo]->s[state].pi = pi;

      /* fixed is a propety of the distribution and optional */
      fixed = getIntAttribute(elem, (const xmlChar *)"fixed", &error);
      if (error)
	fixed = 0;
      else
	f->model.d[modelNo]->s[state].fix = fixed;

      /* order is optional for discrete */
      if (f->modelType & GHMM_kHigherOrderEmissions) {
	order = getIntAttribute(elem, (const xmlChar *)"order", &error);
	if (error)
	  order = 0;
	f->model.d[modelNo]->order[state] = order;
      } else
	order = 0;

      s = (char *)xmlNodeGetContent(elem);
      ARRAY_MALLOC(emissions, pow(f->model.d[modelNo]->M, order+1));
      parseCSVList(s, pow(f->model.d[modelNo]->M, order+1), emissions);
      f->model.d[modelNo]->s[state].b = emissions;
      m_free(s);
    }

    /* ======== continuous state ========================================== */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"mixture"))) {
      assert(f->modelType & GHMM_kContinuousHMM);
      M = 0;
      child = elem->children;
      while (child != NULL) {
        if ((!xmlStrcmp(child->name, (const xmlChar *)"normal")) || 
            (!xmlStrcmp(child->name, (const xmlChar *)"normalTruncatedLeft")) ||
            (!xmlStrcmp(child->name, (const xmlChar *)"normalTruncatedRight")) ||
            (!xmlStrcmp(child->name, (const xmlChar *)"uniform"))){
          M ++;
          
        }
        child = child->next;
      }
      ghmm_c_state_alloc(f->model.c[modelNo]->s + state, M, inDegree[state], outDegree[state], f->model.c[modelNo]->cos);

      f->model.c[modelNo]->M = M;
           
      fixed = getIntAttribute(elem, (const xmlChar *)"fixed", &error);
      if (error) /* optional atribute not defined */          
        fixed = 0;
      
      f->model.c[modelNo]->s[state].fix = fixed;
      f->model.c[modelNo]->s[state].M = M;
      f->model.c[modelNo]->s[state].pi = pi;
      
      child = elem->children;
      
      i = 0;
      while (child != NULL) {  
        if ((!xmlStrcmp(child->name, (const xmlChar *)"normal"))) {
          f->model.c[modelNo]->s[state].mue[i] = getDoubleAttribute(child, (const xmlChar *)"mean", &error);
          f->model.c[modelNo]->s[state].u[i] = getDoubleAttribute(child, (const xmlChar *)"variance", &error);
          f->model.c[modelNo]->s[state].density[i] = (ghmm_density_t)normal;
          aprox = getIntAttribute(child, (const xmlChar *)"aprox", &error);
          if (aprox){
	    f->model.c[modelNo]->s[state].density[i] = (ghmm_density_t)normal_approx;
          }else{
            f->model.c[modelNo]->s[state].density[i] = (ghmm_density_t)normal;
          }          
          fixed = getIntAttribute(child, (const xmlChar *)"fixed", &error);
          if (error)
            fixed = 0;
          f->model.c[modelNo]->s[state].mixture_fix[i] = fixed;
          prior = getDoubleAttribute(child, (const xmlChar *)"prior", &error);
          if (error)
            prior = 1.0;
	  f->model.c[modelNo]->s[state].c[i] = prior;
          i++;                  
        }
        if ((!xmlStrcmp(child->name, (const xmlChar *)"normalTruncatedLeft"))) {
          f->model.c[modelNo]->s[state].mue[i] = getDoubleAttribute(child, (const xmlChar *)"mean", &error);
          f->model.c[modelNo]->s[state].u[i] = getDoubleAttribute(child, (const xmlChar *)"variance", &error);
          f->model.c[modelNo]->s[state].a[i] = getDoubleAttribute(child, (const xmlChar *)"min", &error);
          f->model.c[modelNo]->s[state].density[i] = (ghmm_density_t)normal_left;
          fixed = getIntAttribute(child, (const xmlChar *)"fixed", &error);
          if (error)
            fixed = 0;
          f->model.c[modelNo]->s[state].mixture_fix[i] = fixed;
          prior = getDoubleAttribute(child, (const xmlChar *)"prior", &error);
          if (error)
            prior = 1.0;
	  f->model.c[modelNo]->s[state].c[i] = prior;
          i++;          
        }
        if ((!xmlStrcmp(child->name, (const xmlChar *)"normalTruncatedRight"))) {
          f->model.c[modelNo]->s[state].mue[i] = getDoubleAttribute(child, (const xmlChar *)"mean", &error);
          f->model.c[modelNo]->s[state].u[i] = getDoubleAttribute(child, (const xmlChar *)"variance", &error);
          f->model.c[modelNo]->s[state].a[i] = getDoubleAttribute(child, (const xmlChar *)"max", &error);
          f->model.c[modelNo]->s[state].density[i] = (ghmm_density_t)normal_right;
          fixed = getIntAttribute(child, (const xmlChar *)"fixed", &error);
          if (error)
            fixed = 0;
          f->model.c[modelNo]->s[state].mixture_fix[i] = fixed;
          prior = getDoubleAttribute(child, (const xmlChar *)"prior", &error);
          if (error)
            prior = 1.0;
	  f->model.c[modelNo]->s[state].c[i] = prior;
          i++;  
        }
        if ((!xmlStrcmp(child->name, (const xmlChar *)"uniform"))) {
          f->model.c[modelNo]->s[state].mue[i] = getDoubleAttribute(child, (const xmlChar *)"max", &error);
          f->model.c[modelNo]->s[state].u[i] = getDoubleAttribute(child, (const xmlChar *)"min", &error);
          f->model.c[modelNo]->s[state].density[i] = (ghmm_density_t)uniform;
          fixed = getIntAttribute(child, (const xmlChar *)"fixed", &error);
          if (error)
            fixed = 0;
          f->model.c[modelNo]->s[state].mixture_fix[i] = fixed;
          prior = getDoubleAttribute(child, (const xmlChar *)"prior", &error);
          if (error)
            prior = 1.0;
	  f->model.c[modelNo]->s[state].c[i] = prior;
          i++;      
        }

        child = child->next;         
      }
    }

    /* ======== pair hmm state ============================================ */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"pair"))) {
    }

    /* -------- background name  ------------------------------------------ */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"backgroundName"))) {
      
      assert(f->modelType & GHMM_kBackgroundDistributions);

      s = (char *)xmlNodeGetContent(elem);

      for (i=0; i<f->model.d[modelNo]->bp->n; i++) {
	if (!strcmp(s, f->model.d[modelNo]->bp->name[i])) {
	  if (order != f->model.d[modelNo]->bp->order[i]) {
	    estr = ighmm_mprintf(NULL, 0, "order of background %s and state %d does not match",
				 f->model.d[modelNo]->bp->name[i], state);
	    GHMM_LOG(LERROR, estr);
	    m_free(estr);
	    goto STOP;
	  } else {
	    f->model.d[modelNo]->background_id[state] = i;
	    break;
	  }
	}
      }
      if (i == f->model.d[modelNo]->bp->n) {
	estr = ighmm_mprintf(NULL, 0, "can't find background with name %s in state %d",
			     s, state);
	GHMM_LOG(LERROR, estr);
	m_free(estr);
	goto STOP;
      }
      m_free(s);
    }

    /* -------- tied to --------------------------------------------------- */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"tiedTo"))) {

      assert(f->modelType & GHMM_kTiedEmissions);

      s = (char *)xmlNodeGetContent(elem);
      tied = strtoul(s, serror, 10);
      if (s != *serror && !*serror && state>=tied) {
	f->model.d[modelNo]->tied_to[state] = tied;
	if (f->model.d[modelNo]->tied_to[tied] != tied) {
	  estr = ighmm_mprintf(NULL, 0, "state %d not tied to tie group leader", state);
	  GHMM_LOG(LERROR, estr);
	  m_free(estr);
	  goto STOP;
	}
      } else {
	estr = ighmm_mprintf(NULL, 0, "state %d tiedTo is invalid", state);
	GHMM_LOG(LERROR, estr);
	m_free(estr);
	goto STOP;
      }
      m_free(s);
    }

    /* -------- position for graphical editing ---------------------------- */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"position"))) {
      switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
      case GHMM_kDiscreteHMM:
        f->model.d[modelNo]->s[state].xPosition = getIntAttribute(cur, (const xmlChar *)"x", &error);
        f->model.d[modelNo]->s[state].yPosition = getIntAttribute(cur, (const xmlChar *)"y", &error);
	break;
      case (GHMM_kDiscreteHMM + GHMM_kPairHMM):
        f->model.dp[modelNo]->s[state].xPosition = getIntAttribute(cur, (const xmlChar *)"x", &error);
        f->model.dp[modelNo]->s[state].yPosition = getIntAttribute(cur, (const xmlChar *)"y", &error);
	break;
      case GHMM_kContinuousHMM:
        f->model.c[modelNo]->s[state].xPosition = getIntAttribute(cur, (const xmlChar *)"x", &error);
        f->model.c[modelNo]->s[state].yPosition = getIntAttribute(cur, (const xmlChar *)"y", &error);
	break;
      default:
	GHMM_LOG(LCRITIC, "invalid modelType");} 
    }

    elem = elem->next;
  }

  return 0;
STOP:
  m_free(s);
  m_free(desc);
  m_free(emissions)
  return -1;
#undef CUR_PROC
}

/*===========================================================================*/
static int parseSingleTransition(xmlDocPtr doc, xmlNodePtr cur, fileData_s * f, int modelNo) {
#define CUR_PROC "parseTransition"

  int retval=-1;
  int source, target, error;
  int in_state, out_state;
  double p;
  char * s;
  xmlNodePtr elem;

  assert((f->modelType & GHMM_kTransitionClasses) == 0);

  source = getIntAttribute(cur, (const xmlChar *)"source", &error);
  target = getIntAttribute(cur, (const xmlChar *)"target", &error);

  elem = cur->children;
  while (elem!=NULL) {
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"probability"))) {
      s = (char *)xmlNodeGetContent(elem);
      p = atof(s);
      m_free(s)
      break;
    }
    elem = elem->next;
  }

  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
  case GHMM_kDiscreteHMM:
    out_state = f->model.d[modelNo]->s[source].out_states++;
    in_state  = f->model.d[modelNo]->s[target].in_states++;
    f->model.d[modelNo]->s[source].out_id[out_state] = target;
    f->model.d[modelNo]->s[source].out_a[out_state]  = p;
    f->model.d[modelNo]->s[target].in_id[in_state]   = source;
    f->model.d[modelNo]->s[target].in_a[in_state]    = p;
    break;
  case (GHMM_kDiscreteHMM + GHMM_kPairHMM):
    out_state = f->model.dp[modelNo]->s[source].out_states++;
    in_state  = f->model.dp[modelNo]->s[target].in_states++;
    f->model.dp[modelNo]->s[source].out_id[out_state]   = target;
    f->model.dp[modelNo]->s[source].out_a[out_state][0] = p;
    f->model.dp[modelNo]->s[target].in_id[in_state]     = source;
    f->model.dp[modelNo]->s[target].in_a[in_state][0]   = p;
    break;
  case GHMM_kContinuousHMM:
    out_state = f->model.c[modelNo]->s[source].out_states++;
    in_state  = f->model.c[modelNo]->s[target].in_states++;
    f->model.c[modelNo]->s[source].out_id[out_state]   = target;
    f->model.c[modelNo]->s[source].out_a[0][out_state] = p;
    f->model.c[modelNo]->s[target].in_id[in_state]     = source;
    f->model.c[modelNo]->s[target].in_a[0][in_state]   = p;
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");}
  
  retval = 0;

  return retval;
#undef CUR_PROC
}

/*===========================================================================*/
static int parseMultipleTransition(xmlDocPtr doc, xmlNodePtr cur, fileData_s * f, int modelNo) {
#define CUR_PROC "parseTransition"

  int i, retval=-1;
  int source, target, error, nrTransitionClasses;
  int in_state, out_state;
  double * probs;
  char * s;
  xmlNodePtr elem;

  assert(f->modelType & GHMM_kTransitionClasses);

  source = getIntAttribute(cur, (const xmlChar *)"source", &error);
  target = getIntAttribute(cur, (const xmlChar *)"target", &error);

  elem = cur->children;
  while (elem!=NULL) {
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"probability"))) {
     
      s = (char *)xmlNodeGetContent(elem);
      ARRAY_MALLOC(probs, nrTransitionClasses);
      parseCSVList(s, nrTransitionClasses, probs);
      m_free(s);
      break;
    }
    elem = elem->next;
  }

  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
		+ GHMM_kPairHMM + GHMM_kContinuousHMM)) {
  case (GHMM_kDiscreteHMM + GHMM_kTransitionClasses):
    out_state = f->model.d[modelNo]->s[source].out_states++;
    in_state  = f->model.d[modelNo]->s[target].in_states++;
    f->model.d[modelNo]->s[source].out_id[out_state] = target;
/*XXX    f->model.d[modelNo]->s[source].out_a[out_state]  = probs; */
    f->model.d[modelNo]->s[target].in_id[in_state]   = source;
/*XXX     f->model.d[modelNo]->s[target].in_a[in_state]    = probs; */
    break;
  case (GHMM_kDiscreteHMM + GHMM_kPairHMM + GHMM_kTransitionClasses):
    out_state = f->model.dp[modelNo]->s[source].out_states++;
    in_state  = f->model.dp[modelNo]->s[target].in_states++;
    f->model.dp[modelNo]->s[source].out_id[out_state] = target;
    f->model.dp[modelNo]->s[source].out_a[out_state]  = probs;
    f->model.dp[modelNo]->s[target].in_id[in_state]   = source;
    f->model.dp[modelNo]->s[target].in_a[in_state]    = probs;
    break;
  case (GHMM_kContinuousHMM + GHMM_kTransitionClasses):
    out_state = f->model.c[modelNo]->s[source].out_states++;
    in_state  = f->model.c[modelNo]->s[target].in_states++;
    f->model.c[modelNo]->s[source].out_id[out_state] = target;
    f->model.c[modelNo]->s[target].in_id[in_state]   = source;
    for (i=0; i<nrTransitionClasses; i++) {
      f->model.c[modelNo]->s[source].out_a[i][out_state] = probs[i];
      f->model.c[modelNo]->s[target].in_a[i][in_state]   = probs[i];
    }
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");}

  retval = 0;

STOP:
  m_free(probs);
  return retval;
#undef CUR_PROC
}


/*===========================================================================*/
static int parseHMM(fileData_s * f, xmlDocPtr doc, xmlNodePtr cur, int modelNo) {
#define CUR_PROC "parseHMM"
  char * estr;

  xmlNodePtr child;

  int i, id, error;
  int source, target;

  int N = 0;
  int nrBackgrounds, M=-1;
  int * inDegree = NULL;  
  int * outDegree = NULL;
  int cos;
  float prior;

  int modeltype=0;
  char * mt;
  char * modelname;

  int * bg_orders = NULL;
  double * * bg_ptr = NULL;

  alphabet_s * alfa, * * alphabets;
  int nrAlphabets=0;

  child = cur->children;

   /*xmlElemDump(stdout, doc, cur);*/

  /* parse HMM for counting */
  GHMM_LOG(LINFO, "parseHMM to count ");
  while (child != NULL) {
    /* ========== ALPHABETS ================================================ */
    if ((!xmlStrcmp(child->name, (const xmlChar *)"alphabet"))) {
      if (alphabets)
	ARRAY_MALLOC(alphabets, MAX_ALPHABETS);

      alfa = parseAlphabet(doc, child, f);
      if (alfa && nrAlphabets<MAX_ALPHABETS) {
	alphabets[nrAlphabets++] = alfa;
      } else {
	GHMM_LOG(LERROR, "Error in parsing alphabets.");
	goto STOP;
      }
    }

    /* ========== NODES ==================================================  */
    if ((!xmlStrcmp(child->name, (const xmlChar *)"state"))) {
      id = getIntAttribute(child, (const xmlChar *)"id", &error);
      if (error || id!=N) {
        GHMM_LOG(LERROR, "non consecutive node ids");
        goto STOP;
      }
      N++;
    }

    /* ========== EDGES ==================================================  */
    if ((!xmlStrcmp(child->name, (const xmlChar *)"transition"))) {
      
      if (inDegree == NULL) {
        ARRAY_CALLOC(inDegree, N);
        ARRAY_CALLOC(outDegree, N);
      }

      source = getIntAttribute(child, (const xmlChar *)"source", &error);
      if (error || source<0 || source>N) {
	estr = ighmm_mprintf(NULL, 0, "source (%d) node not existing (%d)", source, error);
	GHMM_LOG(LERROR, estr);
	m_free(estr);
        goto STOP;
      }

      target = getIntAttribute(child, (const xmlChar *)"target", &error);
      if (error || target<0 || target>N) {
	estr = ighmm_mprintf(NULL, 0, "target (%d) node not existing (%d)", target, error);
        GHMM_LOG(LERROR, estr);
	m_free(estr);
        goto STOP;
      }

      inDegree[target]++;
      outDegree[source]++;
    }
    /* ========== BACKGROUND DISTRIBUTIONS ================================  */
    if ((!xmlStrcmp(child->name, (const xmlChar *)"background")))
      nrBackgrounds++;

    child = child->next;
  }

  estr = ighmm_mprintf(NULL, 0, "Found HMM with %d states\n", N);
  GHMM_LOG(LDEBUG, estr);
  m_free(estr);
  for (i=0; i<N; i++) {
    estr = ighmm_mprintf(NULL, 0, "  %d\t%d\n", inDegree[i], outDegree[i]);
    GHMM_LOG(LDEBUG, estr);
    m_free(estr);
  }

  /* starting real parsing */
  modelname = (char *)getXMLCharAttribute(cur, (const xmlChar *)"name", &error);
  mt = (char *)getXMLCharAttribute(cur, (const xmlChar *)"type", &error);
  modeltype = parseModelType(mt, strlen(mt));
  f->modelType = modeltype;

  /* if firtst model, initialize model structures */
  if ( modelNo == 0){
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
    case GHMM_kDiscreteHMM:
      ARRAY_CALLOC(f->model.d,f->noModels);
      break;
    case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
      ARRAY_CALLOC(f->model.ds,f->noModels);
      break;
    case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
    case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
      ARRAY_CALLOC(f->model.dp,f->noModels);  
      break;
    case GHMM_kContinuousHMM:
    case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):    
      ARRAY_CALLOC(f->model.c,f->noModels);  
      break;
    break;
    default:
      GHMM_LOG(LCRITIC, "invalid modelType");
    }       
  }

  /* allocating the different models */
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
  case GHMM_kDiscreteHMM:
    assert(nrAlphabets == 1);
    M = alphabets[0]->size;
    f->model.d[modelNo] = ghmm_dmodel_calloc(M, N, modeltype, inDegree, outDegree);
    f->model.d[modelNo]->alphabet = alphabets[0];
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    f->model.d[modelNo] = NULL;
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    /* f->model.dp[modelNo] = NULL; XXX*/
    break;
  case GHMM_kContinuousHMM:
  case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):    
    f->model.c[modelNo] = ghmm_cmodel_calloc(N,modeltype);
    prior = getDoubleAttribute(cur, (const xmlChar *)"prior", &error);
    if (error) /* optional atribute not defined */          
        prior = 0.0;      
    f->model.c[modelNo]->prior = prior;
    cos = getIntAttribute(cur, (const xmlChar *)"transitionClasses", &error);
    if (error) /* optional atribute not defined */          
        cos = 1;      
    f->model.c[modelNo]->cos = cos;
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");
  }

  /* allocating background distributions for approtiate models */
  if (modeltype & GHMM_kBackgroundDistributions) { 
    switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			    + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
    case GHMM_kDiscreteHMM:
      ARRAY_CALLOC(bg_orders, N);
      ARRAY_CALLOC(bg_ptr, N);
      f->model.d[modelNo]->bp = ghmm_d_background_alloc(nrBackgrounds, M, bg_orders, bg_ptr);
      ARRAY_CALLOC(f->model.d[modelNo]->bp->name, N);
      f->model.d[modelNo]->bp->n = 0;
      break;
    default:
      GHMM_LOG(LERROR, "invalid modelType");}
  }

  child = cur->xmlChildrenNode;

  /* parse HMM for real */
  while (child != NULL) {

    /* ========== LABEL ALPHABETS ========================================== */
    if ((!xmlStrcmp(child->name, (const xmlChar *)"labelAlphabet"))) {
      alfa = parseAlphabet(doc, child, f);
      if (alfa) {
	f->model.d[modelNo]->labelAlphabet = alfa;
      } else {
	GHMM_LOG(LERROR, "Error in parsing alphabets.");
	goto STOP;
      }
    }
    
    if ((!xmlStrcmp(child->name, (const xmlChar *)"background"))) {
      if (modeltype & GHMM_kBackgroundDistributions) {
	parseBackground(doc, child, f, modelNo);
      } else {
	GHMM_LOG(LWARN, "Ignoring background distribution.");
      }
    }
    if ((!xmlStrcmp(child->name, (const xmlChar *)"state"))) {
      parseState(doc, child, f, inDegree, outDegree, modelNo);
    }
    if ((!xmlStrcmp(child->name, (const xmlChar *)"transition"))) {
      if (modeltype & GHMM_kTransitionClasses)
	parseMultipleTransition(doc, child, f, modelNo);
      else
	parseSingleTransition(doc, child, f, modelNo);
    }
    child = child->next;
  }

  /* freeing temporary data */
  m_free(inDegree);
  m_free(outDegree);  
  return 0;  
STOP:
  free(inDegree);
  free(outDegree);
  free(bg_orders);
  free(bg_ptr);
  free(alphabets);
  free(f);
  return -1;
#undef CUR_PROC
}


/*===========================================================================*/
fileData_s * parseHMMDocument(const char *filename) {
#define CUR_PROC "parseHMMDocument"

  xmlParserCtxtPtr ctxt; /* the parser context */
  xmlDocPtr doc; /* the resulting document tree */
  xmlNodePtr cur, child;
  int modelNo = 0;
  int error;

  char * str;
  fileData_s * filedata = NULL;

  /* create a parser context */
  ctxt = xmlNewParserCtxt();
  if (ctxt == NULL) {
    GHMM_LOG(LERROR, "Failed to allocate parser context");
    return NULL;
  }
  /* parse the file, activating the DTD validation option */
  doc = xmlCtxtReadFile(ctxt, filename, NULL, XML_PARSE_DTDVALID);
  /* check if parsing suceeded */
  if (doc == NULL) {
    str = ighmm_mprintf(NULL, 0, "Failed to parse %s", filename);
    GHMM_LOG(LERROR, str);
    m_free(str);
  } else {
    /* check if validation suceeded */
    if (ctxt->valid == 0) {
      str = ighmm_mprintf(NULL, 0, "Failed to validate %s", filename);
      GHMM_LOG(LERROR, str);
      m_free(str);
    } else {
      /* checking the root node, creating the file structure and iteration over all HMMs */
      cur = xmlDocGetRootElement(doc);
      if ((!xmlStrcmp(cur->name, (const xmlChar *)"mixture"))) {
        ARRAY_CALLOC(filedata, 1);        
        filedata->noModels = getIntAttribute(cur, (const xmlChar *)"noComponents", &error);                
        child = cur->children;
        while (child!=NULL) {
          if ((!xmlStrcmp(child->name, (const xmlChar *)"HMM"))) {
            if (modelNo > filedata->noModels){
              str = ighmm_mprintf(NULL, 0, "The mixture has more models than defined");
	      GHMM_LOG(LERROR, str);
	      m_free(str);
            }else{
              parseHMM(filedata,doc, child, modelNo);
              modelNo++;
            }
          }
          child=child->next; 
        }
        if (modelNo < filedata->noModels){
          str = ighmm_mprintf(NULL, 0, "The mixture has less models than defined");
	  GHMM_LOG(LERROR, str);
	  m_free(str);
        }
      }else{
        str = ighmm_mprintf(NULL, 0, "The file does not contains the appropriate root %s", filename);
	GHMM_LOG(LERROR, str);
	m_free(str);
      }
    }

    /* free up the resulting document */
    xmlFreeDoc(doc);
  }
  /* free up the parser context */
  xmlFreeParserCtxt(ctxt);

  return filedata;
STOP:
  /*do this */ 
  return NULL; 
#undef CUR_PROC
}
