/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/ghmm/xmlio.h
  created      : DATE: May 2001
  $Id$

__copyright__

*******************************************************************************/

#ifndef XMLIO_H
#define XMLIO_H
#include <expat.h>

/*
  identify different handler_types
 */
typedef enum {
  document_handler,
  vertex_handler,
  default_handler
} object_handler_type;


/* these object handlers are stored in a double linked list
   the first element is an document handler, that can handle all events.

   - when a new element in the document is found a special handler is created and appended to the list.
   - if an event can't be handled by the actual object handler, it is propagated through the list.
   - recieved data are stored by the actual handler or the parents.
   - if an element ends the element user_data are given to the parent handler and the handler is deleted

*/
typedef struct object_handler object_handler;

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
  
  /* parser */
  XML_Parser parser;

  /* its data */
  void* handler_data;

  /* its functions */
  object_handler_functions functions;
};

object_handler*
default_create_handler_object(XML_Parser parser,
			      const XML_Char *name,
			      const XML_Char **atts);

void
default_delete_handler_object(object_handler* handler);

void
default_set_handler(object_handler* handler);

object_handler*
default_push_handler(object_handler* old_handler,
		     object_handler* new_handler);

object_handler*
default_pop_handler(object_handler* handler);

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

int
ghmm_xml_parse(const char* filename);

#endif /* XMLIO_H */



