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

typedef void (*ChildDataHandler)(void* ChildData, object_handler_type ChildType);

typedef struct 
{
  XML_StartElementHandler StartHandler;
  XML_EndElementHandler EndHandler;
  XML_CharacterDataHandler CharacterHandler;
  ChildDataHandler ChildDataReciever;
} object_handler_functions;

/* these object handlers are stored in a double linked list
   the first element is an document handler, that can handle all events.

   - when a new element in the document is found a special handler is created and appended to the list.
   - if an event can't be handled by the actual object handler, it is propagated through the list.
   - recieved data are stored by the actual handler or the parents.
   - if an element ends the element user_data are given to the parent handler and the handler is deleted

*/
typedef struct object_handler object_handler;

/* */
struct object_handler
{
  /* Double linked list */
  object_handler *prec;
  object_handler *succ;

  /* type of this handler */
  object_handler_type type;
  
  /* its data */
  void* handler_data;

  /* its functions */
  object_handler_functions functions;
};

#endif /* XMLIO_H */
