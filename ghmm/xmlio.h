/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/ghmm/xmlio.h
  created      : DATE: May 2001
  $Id$

__copyright__

*******************************************************************************/

#ifndef XMLIO_H
#define XMLIO_H
#include <stdio.h>
#include <expat.h>

/*
  identify different handler_types
 */
typedef enum {
  document_object,
  vertex_object,
  default_object
} object_handler_type;


/* these object handlers are stored in a double linked list
   the first element is an document handler, that can handle all events.

   - when a new element in the document is found a special handler is created and appended to the list.
   - if an event can't be handled by the actual object handler, it is propagated through the list.
   - recieved data are stored by the actual handler or the parents.
   - if an element ends the element user_data are given to the parent handler and the handler is deleted

*/
typedef struct object_handler object_handler;
typedef struct document_handler document_handler;

/*
  evaluates contents of parsed elements
 */
typedef void (*ChildDataHandler)(object_handler* ChildData);

typedef struct 
{
  XML_StartElementHandler StartHandler;
  XML_EndElementHandler EndHandler;
  XML_CharacterDataHandler CharacterHandler;
  ChildDataHandler ChildDataReciever;
} object_handler_functions;


/* */
struct object_handler
{
  /* Double linked list */
  object_handler *prec;
  object_handler *succ;

  /* type of this handler */
  object_handler_type type;
  
  /* document information */
  document_handler* dh;

  /* its data */
  void* handler_data;

  /* its functions */
  object_handler_functions functions;
};

object_handler*
default_create_handler_object(document_handler* dh,
			      const XML_Char *name,
			      const XML_Char **atts);

void
default_delete_handler_object(object_handler* handler);

void
default_StartElement_handler(void *userData,
			     const XML_Char *name,
			     const XML_Char **atts);

void
default_CharacterData_handler(void *userData,
			      const XML_Char *s,
			      int len);

void
default_EndElement_handler(void *userData,
			   const XML_Char *name);

void
default_ChildData_handler(object_handler* handler);

/*
  more abstract
*/

object_handler*
create_object_handler(document_handler* dh,
		      const object_handler_type type,
		      const XML_Char *name,
		      const XML_Char **atts);

void
delete_object_handler(object_handler* handler);

struct  document_handler{
  FILE* xml_file;
  XML_Parser* parser;
  object_handler root;
  object_handler* last;
};


object_handler*
document_handler_push_object_handler(document_handler* dh,
				     object_handler* new_handler);


object_handler*
document_handler_pop_object_handler(document_handler* dh);

document_handler*
create_document_handler(const char* filename);

void
document_handler_set_handler(object_handler* handler);

int
parse_document(document_handler* dh);

int
ghmm_xml_parse(const char* filename);

#endif /* XMLIO_H */









