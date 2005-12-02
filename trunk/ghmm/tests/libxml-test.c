#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include <ghmm/ghmm_internals.h>

int getIntAttribute(xmlNodePtr node, const xmlChar *name, int *error)
{
  xmlChar *attr;
  int value = 0;

  if( (attr = xmlGetProp(node, name)) != NULL) {
    value = atoi((char*)attr);
    xmlFree(attr);
    *error = 0;
  } else {
    *error = 1;
  }
  return value;
}

/* Caller owns return value */
xmlChar* getXMLCharAttribute(xmlNodePtr node, const xmlChar *name, int *error)
{
  xmlChar *attr;

  if( (attr = xmlGetProp(node, name)) != NULL) {
    *error = 0;
    return attr;
  } else {
    *error = 1;
    return NULL;
  }
}

static void 
parseNode(xmlDocPtr doc, xmlNodePtr cur)
{
  int error;
  xmlChar *stuff;
  xmlNode *data = NULL;
  xmlChar* dataKey;
  char *s, *tok;

  /* Handle Node Data */
  for( data = cur->xmlChildrenNode; data != NULL; data = data->next) {

    if(data->type == XML_ELEMENT_NODE) {
      if( !xmlStrcmp(data->name, (const xmlChar*)"data")) {
	dataKey = getXMLCharAttribute(data,(const xmlChar*)"key",&error);
	printf("  Data key=%s\n", (char*)dataKey);
	    
	/* Collect the value */
	stuff = xmlNodeListGetString(doc, data->xmlChildrenNode, 1);

	if (strcmp( (char*)dataKey, "emissions" ) == 0) {
	  /* parse emission vector */
	  s = (char *) stuff;
	  tok = strtok(s, ",");
	  do {
	    printf( "  %s ", tok);
	    tok = strtok(NULL, ",");
	  } while ( tok != NULL );
	  printf( "\n" );
	} else {
	  printf("  value=%s\n", (char*)stuff);
	}

	xmlFree(stuff);
	xmlFree(dataKey);
      }
    }
  }
}

static void
parseGraph(xmlDocPtr doc, xmlNodePtr cur)
{
#define CUR_PROC "parseGraph"

  xmlChar *stuff;
  xmlNode *data = NULL;
  xmlChar* dataKey;
  int i, id, error;
  int source, target;

  int N = 0; /* We don't follow C-style indexing strategy */ 
  int *inDegree = NULL;  
  int *outDegree = NULL;

  cur = cur->xmlChildrenNode;

  GHMM_LOG(LINFO, "parseGraph");
  while( cur != NULL) {

    /* ========== NODES ==================================================  */
    if( (!xmlStrcmp(cur->name, (const xmlChar*)"node"))) {      
      id = getIntAttribute(cur, (const xmlChar*)"id", &error); /*\ufffdGet the node ID */
      if(!error)
	printf("node id = %d\n", id);

      N=id+1;
      parseNode(doc, cur);
    }
    
    /* ========== EDGES ==================================================  */
    if( (!xmlStrcmp(cur->name, (const xmlChar*)"edge"))) {      

      if (inDegree == NULL) {
	inDegree  = calloc(N+1, sizeof(int));
	outDegree = calloc(N+1, sizeof(int));
      }

      source = getIntAttribute(cur, (const xmlChar*)"source", &error); /*\ufffdGet the node ID */
      if(!error)
	printf("source = %d\n", source);
      else
	return;

      target = getIntAttribute(cur, (const xmlChar*)"target", &error); /*\ufffdGet the node ID */
      if(!error)
	printf("target = %d\n", target);      
      else
	return;
      
      inDegree[target] += 1;
      outDegree[source] += 1;
    }
    cur = cur->next;
  }

  printf("Found HMM with %d states\n", N);
  for(i = 1; i <= N; i++) {
    printf("  %d\t%d\n", inDegree[i], outDegree[i]);
  }
  
  free(inDegree);
  free(outDegree);
#undef CUR_PROC
}


static void
parseHMMDocument(char *docname) {
  
  xmlDocPtr doc;
  xmlNodePtr cur;

  doc = xmlParseFile(docname);

  if(doc == NULL) {
    fprintf(stderr, "Document %s not parsed successfully.\n", docname);
    return;
  }
  
  cur = xmlDocGetRootElement(doc);

  if(cur == NULL) {
    fprintf(stderr, "Document %s is empty\n", docname); 
    xmlFreeDoc(doc);
    return;
  }

  if( xmlStrcmp(cur -> name, (const xmlChar*)"HMM")) {
    fprintf(stderr, "Document %s is not of type 'HMM'\n", docname); 
    xmlFreeDoc(doc);
    return;
  }
  else {
    parseGraph(doc, cur);
  }

  xmlFreeDoc(doc);
  return;
}


int
main(int argc, char **argv) {

  char *docname;

  if(argc <= 1) {
    printf("Usage: %s docname.xml", argv[0]);
    return(0);
  }
  docname = argv[1];
  parseHMMDocument(docname);
  
  return(1);
}
