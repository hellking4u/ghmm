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

#include <ghmm/ghmm.h>
#include <ghmm/ghmm_internals.h>
#include <ghmm/mes.h>
#include <ghmm/mprintf.h>
#include <ghmm/model.h>
#include <ghmm/smodel.h>
#include <ghmm/pmodel.h>
#include <ghmm/sdmodel.h>

/* we should not need more */
#define MAX_ALPHABETS 2 

struct alphabet_s {
  int id;
  char * description;
  unsigned int size;
  char * * symbols;
};
typedef struct alphabet_s alphabet_s;

struct fileData_s {

  int modelType;

  union {
    ghmm_cmodel * c;
    ghmm_dmodel * d;
    ghmm_dpmodel * dp;
    ghmm_dsmodel * ds;
  } model;

  ghmm_cmodel * aux;

  unsigned int nrAlphabets;
  alphabet_s * * alphabets;
  alphabet_s * labelAlphabet;

  int * xPosition;
  int * yPosition;

  

};
typedef struct fileData_s fileData_s;


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
static int parseBackground(xmlDocPtr doc, xmlNodePtr cur, fileData_s * f) {
#define CUR_PROC "parseBackground"

  int error, order;
  int bgNr;
  double * b;
  char * s;

  assert(f->modelType & GHMM_kDiscreteHMM);

  bgNr = f->model.d->bp->n++;

  /* get order */
  order = getIntAttribute(cur, (const xmlChar *)"order", &error);
  if (error)
    order=0;
  else if (order && !(f->modelType & GHMM_kHigherOrderEmissions)) {
    GHMM_LOG(LERROR, "background distribution has order > 0, but model is not higher order");
    goto STOP;
  }
  f->model.d->bp->order[bgNr] = order;

  /* get name */
  s = (char *)getXMLCharAttribute(cur, (const xmlChar *)"name", &error);
  f->model.d->bp->name[bgNr] = s;

  /* get distribution */
  s = (char *)xmlNodeGetContent(cur);

  ARRAY_MALLOC(b, pow(f->model.d->bp->m, order+1));
  if (-1 !=  parseCSVList(s, pow(f->model.d->bp->m, order+1), b))
    f->model.d->bp->b[bgNr] = b;
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
static int parseState(xmlDocPtr doc, xmlNodePtr cur, fileData_s * f, int * inDegree, int * outDegree) {
#define CUR_PROC "parseState"

  int i, error, order=-28, state=-1442, fixed=-985, tied=-9354, M, cos, aprox;
  double pi;
  double * emissions;
  char * desc, * s, * estr, * * serror;


  xmlNodePtr elem, child;

  state = getIntAttribute(cur, (const xmlChar *)"id", &error);
  pi = getDoubleAttribute(cur, (const xmlChar *)"inital", &error);
  desc = (char *)getXMLCharAttribute(cur, (const xmlChar *)"desc", &error);

  elem = cur->children;
  while (elem!=NULL) {
    /* ======== silent state ============================================== */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"silent"))) {
      if (f->modelType & GHMM_kDiscreteHMM)
	f->model.d->silent[state] = 1;
    }

    /* ======== discrete state (possible higher order) ==================== */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"discrete"))) {
      assert((f->modelType & GHMM_kDiscreteHMM) && ((f->modelType & GHMM_kPairHMM) == 0));

      /* fixed is a propety of the distribution */
      fixed = getIntAttribute(elem, (const xmlChar *)"fixed", &error);
      f->model.d->s[state].fix = fixed;

      /* order is optional for discrete */
      if (f->modelType & GHMM_kHigherOrderEmissions) {
	order = getIntAttribute(elem, (const xmlChar *)"order", &error);
	if (error)
	  order = 0;
	f->model.d->order[state] = order;
      } else
	order = 0;

      s = (char *)xmlNodeGetContent(elem);
      ARRAY_MALLOC(emissions, pow(f->model.d->M, order+1));
      parseCSVList(s, pow(f->model.d->M, order+1), emissions);
      f->model.d->s[state].b = emissions;
      m_free(s);
    }

    /* ======== continuous state ========================================== */
    if ((!xmlStrcmp(elem->name, (const xmlChar *)"mixture"))) {
      assert(f->modelType & GHMM_kContinuousHMM);
      M = 0;
      cos = 1; /*read this*/
      child = elem->children;
      while (child != NULL) {
        M ++;
        child = child->next;
      }
      /*ghmm_c_print(stdout,&f->model.c);*/
      ghmm_c_state_alloc( f->model.c->s + state, M, inDegree[state], outDegree[state], cos);
           

      fixed = getIntAttribute(elem, (const xmlChar *)"fixed", &error);
      if (error) /* optional atribute not defined */          
        fixed = 0;
      
      f->model.d->s[state].fix = fixed;
      f->model.c->s[state].M = M;
      f->model.c->s[state].pi =pi;
      
      child = elem->children;
      
      for(i=0;i<M;i++){  
        if ((!xmlStrcmp(child->name, (const xmlChar *)"normal"))) {
          f->model.c->s[state].mue[i] = getDoubleAttribute(child, (const xmlChar *)"mean", &error);
          f->model.c->s[state].u[i] = getDoubleAttribute(child, (const xmlChar *)"variance", &error);
          f->model.c->s[state].density[i] = (ghmm_density_t)normal;
          aprox = getIntAttribute(child, (const xmlChar *)"aprox", &error);
          if (aprox){
	    f->model.c->s[state].density[i] = (ghmm_density_t)normal_approx;
          }else{
            f->model.c->s[state].density[i] = (ghmm_density_t)normal;
          }                   
        }
        if ((!xmlStrcmp(child->name, (const xmlChar *)"normalTruncatedLeft"))) {
          f->model.c->s[state].mue[i] = getDoubleAttribute(child, (const xmlChar *)"mean", &error);
          f->model.c->s[state].u[i] = getDoubleAttribute(child, (const xmlChar *)"variance", &error);
          f->model.c->s[state].a[i] = getDoubleAttribute(child, (const xmlChar *)"min", &error);
          f->model.c->s[state].density[i] = (ghmm_density_t)normal_left;           
        }
        if ((!xmlStrcmp(child->name, (const xmlChar *)"normalTruncatedRight"))) {
          f->model.c->s[state].mue[i] = getDoubleAttribute(child, (const xmlChar *)"mean", &error);
          f->model.c->s[state].u[i] = getDoubleAttribute(child, (const xmlChar *)"variance", &error);
          f->model.c->s[state].a[i] = getDoubleAttribute(child, (const xmlChar *)"max", &error);
          f->model.c->s[state].density[i] = (ghmm_density_t)normal_right;           
        }
        if ((!xmlStrcmp(child->name, (const xmlChar *)"uniform"))) {
          f->model.c->s[state].mue[i] = getDoubleAttribute(child, (const xmlChar *)"max", &error);
          f->model.c->s[state].u[i] = getDoubleAttribute(child, (const xmlChar *)"min", &error);
          f->model.c->s[state].density[i] = (ghmm_density_t)uniform;             
        }
        fixed = getIntAttribute(child, (const xmlChar *)"fixed", &error);
        f->model.c->s[state].mixture_fix[i] = fixed;
        s = (char *)xmlNodeGetContent(child);
        f->model.c->s[state].c[i] =  atof(s);
        m_free(s);
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

      for (i=0; i<f->model.d->bp->n; i++) {
	if (!strcmp(s, f->model.d->bp->name[i])) {
	  if (order != f->model.d->bp->order[i]) {
	    estr = ighmm_mprintf(NULL, 0, "order of background %s and state %d does not match",
				 f->model.d->bp->name[i], state);
	    GHMM_LOG(LERROR, estr);
	    m_free(estr);
	    goto STOP;
	  } else {
	    f->model.d->background_id[state] = i;
	    break;
	  }
	}
      }
      if (i == f->model.d->bp->n) {
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
	f->model.d->tied_to[state] = tied;
	if (f->model.d->tied_to[tied] != tied) {
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
      f->xPosition[state] = getIntAttribute(cur, (const xmlChar *)"x", &error);
      f->yPosition[state] = getIntAttribute(cur, (const xmlChar *)"y", &error);
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
static int parseSingleTransition(xmlDocPtr doc, xmlNodePtr cur, fileData_s * f) {
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
    out_state = f->model.d->s[source].out_states++;
    in_state  = f->model.d->s[target].in_states++;
    f->model.d->s[source].out_id[out_state] = target;
    f->model.d->s[source].out_a[out_state]  = p;
    f->model.d->s[target].in_id[in_state]   = source;
    f->model.d->s[target].in_a[in_state]    = p;
    break;
  case (GHMM_kDiscreteHMM + GHMM_kPairHMM):
    out_state = f->model.dp->s[source].out_states++;
    in_state  = f->model.dp->s[target].in_states++;
    f->model.dp->s[source].out_id[out_state]   = target;
    f->model.dp->s[source].out_a[out_state][0] = p;
    f->model.dp->s[target].in_id[in_state]     = source;
    f->model.dp->s[target].in_a[in_state][0]   = p;
    break;
  case GHMM_kContinuousHMM:
    out_state = f->model.c->s[source].out_states++;
    in_state  = f->model.c->s[target].in_states++;
    f->model.c->s[source].out_id[out_state]   = target;
    f->model.c->s[source].out_a[0][out_state] = p;
    f->model.c->s[target].in_id[in_state]     = source;
    f->model.c->s[target].in_a[0][in_state]   = p;
    break;
  default:
    GHMM_LOG(LCRITIC, "invalid modelType");}
  
  retval = 0;

  return retval;
#undef CUR_PROC
}

/*===========================================================================*/
static int parseMultipleTransition(xmlDocPtr doc, xmlNodePtr cur, fileData_s * f) {
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
    out_state = f->model.ds->s[source].out_states++;
    in_state  = f->model.ds->s[target].in_states++;
    f->model.ds->s[source].out_id[out_state] = target;
    f->model.ds->s[source].out_a[out_state]  = probs;
    f->model.ds->s[target].in_id[in_state]   = source;
    f->model.ds->s[target].in_a[in_state]    = probs;
    break;
  case (GHMM_kDiscreteHMM + GHMM_kPairHMM + GHMM_kTransitionClasses):
    out_state = f->model.dp->s[source].out_states++;
    in_state  = f->model.dp->s[target].in_states++;
    f->model.dp->s[source].out_id[out_state] = target;
    f->model.dp->s[source].out_a[out_state]  = probs;
    f->model.dp->s[target].in_id[in_state]   = source;
    f->model.dp->s[target].in_a[in_state]    = probs;
    break;
  case (GHMM_kContinuousHMM + GHMM_kTransitionClasses):
    out_state = f->model.ds->s[source].out_states++;
    in_state  = f->model.ds->s[target].in_states++;
    f->model.c->s[source].out_id[out_state] = target;
    f->model.c->s[target].in_id[in_state]   = source;
    for (i=0; i<nrTransitionClasses; i++) {
      f->model.c->s[source].out_a[i][out_state] = probs[i];
      f->model.c->s[target].in_a[i][in_state]   = probs[i];
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
static int parseHMM(xmlDocPtr doc, xmlNodePtr cur) {
#define CUR_PROC "parseHMM"

  char * estr;

  xmlNodePtr child;

  int i, id, error;
  int source, target;

  int N = 0;
  int nrBackgrounds, M=-1;
  int * inDegree = NULL;  
  int * outDegree = NULL;

  int modeltype=0;
  char * mt;
  char * modelname;

  int * bg_orders;
  double * * bg_ptr;

  alphabet_s * alfa;

  fileData_s * f;

  ARRAY_CALLOC(f, 1);

  child = cur->children;

   /*xmlElemDump(stdout, doc, cur);*/

  /* parse HMM for counting */
  GHMM_LOG(LINFO, "parseHMM to count ");
  while (child != NULL) {

    /* ========== ALPHABETS ================================================ */
    if ((!xmlStrcmp(child->name, (const xmlChar *)"alphabet"))) {
      if (!f->alphabets)
	ARRAY_MALLOC(f->alphabets, MAX_ALPHABETS);

      alfa = parseAlphabet(doc, child, f);
      if (alfa && f->nrAlphabets<MAX_ALPHABETS) {
	f->alphabets[f->nrAlphabets++] = alfa;
      } else {
	GHMM_LOG(LERROR, "Error in parsing alphabets.");
	goto STOP;
      }
    }

    /* ========== LABEL ALPHABETS ========================================== */
    if ((!xmlStrcmp(child->name, (const xmlChar *)"labelAlphabet"))) {
      alfa = parseAlphabet(doc, child, f);
      if (alfa) {
	f->labelAlphabet = alfa;
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
  /* allocating space in the file data struct */
  ARRAY_CALLOC(f->xPosition, N);
  ARRAY_CALLOC(f->yPosition, N);

  /* starting real parsing */
  modelname = (char *)getXMLCharAttribute(cur, (const xmlChar *)"name", &error);
  mt = (char *)getXMLCharAttribute(cur, (const xmlChar *)"type", &error);
  modeltype = parseModelType(mt, strlen(mt));
  f->modelType = modeltype;

  /* allocating the different models */
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
  case GHMM_kDiscreteHMM:
    M = f->alphabets[0]->size;
    f->model.d = ghmm_dmodel_calloc(M, N, modeltype, inDegree, outDegree);
    printf("N = %d\n", f->model.d->N);
    break;
  case (GHMM_kDiscreteHMM+GHMM_kTransitionClasses):
    f->model.ds = NULL;
    break;
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM):
  case (GHMM_kDiscreteHMM+GHMM_kPairHMM+GHMM_kTransitionClasses):
    f->model.dp = NULL;
    break;
  case GHMM_kContinuousHMM:
  case (GHMM_kContinuousHMM+GHMM_kTransitionClasses):    
    f->model.c = ghmm_cmodel_calloc(N,modeltype);
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
      f->model.d->bp = ghmm_d_background_alloc(nrBackgrounds, M, bg_orders, bg_ptr);
      ARRAY_CALLOC(f->model.d->bp->name, N);
      f->model.d->bp->n = 0;
      break;
    default:
      GHMM_LOG(LERROR, "invalid modelType");}
  }

  child = cur->xmlChildrenNode;

  /* parse HMM for real */
  while (child != NULL) {

    if ((!xmlStrcmp(child->name, (const xmlChar *)"background"))) {
      if (modeltype & GHMM_kBackgroundDistributions) {
	parseBackground(doc, child, f);
      } else {
	GHMM_LOG(LWARN, "Ignoring background distribution.");
      }
    }
    if ((!xmlStrcmp(child->name, (const xmlChar *)"state"))) {
      parseState(doc, child, f, inDegree, outDegree);
    }
    if ((!xmlStrcmp(child->name, (const xmlChar *)"transition"))) {
      if (modeltype & GHMM_kTransitionClasses)
	parseMultipleTransition(doc, child, f);
      else
	parseSingleTransition(doc, child, f);
    }
    child = child->next;
  }

  /* freeing temporary data */
  m_free(inDegree);
  m_free(outDegree);  
  
  switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			  + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
    case GHMM_kContinuousHMM:
      ghmm_c_print(stdout, f->model.c);
    break;
    case GHMM_kDiscreteHMM:
      ghmm_d_print(stdout, f->model.d);
    default:
      break;
  }

  return 0;
STOP:
  if (inDegree) {
    m_free(inDegree);
    m_free(outDegree);
  }
  m_free(bg_orders);
  m_free(bg_ptr);
  m_free(f->xPosition);
  m_free(f->yPosition);
  m_free(f)
  return -1;
#undef CUR_PROC
}


/*===========================================================================*/
static void parseHMMDocument(const char *filename) {
#define CUR_PROC "parseHMMDocument"

  xmlParserCtxtPtr ctxt; /* the parser context */
  xmlDocPtr doc; /* the resulting document tree */
  xmlNodePtr cur;

  char * str;

  /* create a parser context */
  ctxt = xmlNewParserCtxt();
  if (ctxt == NULL) {
    GHMM_LOG(LCRITIC, "Failed to allocate parser context");
    return;
  }
  /* parse the file, activating the DTD validation option */
  doc = xmlCtxtReadFile(ctxt, filename, NULL, XML_PARSE_DTDVALID);
  /* check if parsing suceeded */
  if (doc == NULL) {
    str = ighmm_mprintf(NULL, 0, "Failed to parse %s", filename);
    GHMM_LOG(LCRITIC, str);
    m_free(str);
  } else {
    /* check if validation suceeded */
    if (ctxt->valid == 0) {
      str = ighmm_mprintf(NULL, 0, "Failed to validate %s", filename);
      GHMM_LOG(LERROR, str);
      m_free(str);
    } else {
      
      cur = xmlDocGetRootElement(doc);

      parseHMM(doc, cur);

    }

    /* free up the resulting document */
    xmlFreeDoc(doc);
  }
  /* free up the parser context */
  xmlFreeParserCtxt(ctxt);
#undef CUR_PROC
}


/*===========================================================================*/
int main(int argc, char **argv) {

  char *docname;

  ghmm_set_loglevel(LDEBUG+1);

  if(argc <= 1) {
    printf("Usage: %s docname.xml", argv[0]);
    return(0);
  }

  /*
   * this initialize the library and check potential ABI mismatches
   * between the version it was compiled for and the actual shared
   * library used.
   */
  LIBXML_TEST_VERSION
    
  docname = argv[1];
  parseHMMDocument(docname);

  /*
   * Cleanup function for the XML library.
   */
  xmlCleanupParser();
  /*
   * this is to debug memory for regression tests
   */
  xmlMemoryDump();
  return(0);
}
