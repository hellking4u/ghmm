/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/ghmm/xmlio.c
  created      : DATE: May 2001
  $Id$

__copyright__

*******************************************************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif /* HAVE_CONFIG_H */

#if defined(__EXPERIMENTAL__) && __EXPERIMENTAL__ == 1
#if defined(HAVE_EXPAT_H) && defined(HAVE_LIBEXPAT)

#include <stdlib.h>
#include <stdio.h>
#include "xmlio.h"

object_handler*
default_create_handler_object(XML_Parser parser,
			      const XML_Char *name,
			      const XML_Char **atts)
{
  /*
    allocate new handler
   */
  object_handler* new_handler;
  new_handler=(object_handler*)malloc(sizeof(object_handler));
  if (new_handler==NULL)
    {
      fprintf(stderr,"can't allocate new object_handler for %s\n",name);
      return NULL;
    }

  /*
    set handler
   */
  new_handler->type=default_handler;
  new_handler->parser=parser;
  new_handler->handler_data=NULL;
  new_handler->functions.StartHandler=&default_StartElement_handler;
  new_handler->functions.EndHandler=&default_EndElement_handler;
  new_handler->functions.CharacterHandler=&default_CharacterData_handler;
  new_handler->functions.ChildDataReciever=&default_ChildData_handler;

  /*
    evaluate name and attributes
   */

  fprintf(stderr,"this is the new %s handler at %p\n",name,new_handler);

  /*
    return this handler
   */
  return new_handler;
}

void
default_delete_handler_object(object_handler* handler)
{
  /*
    do other handler specific things...
   */

  if (handler->handler_data!=NULL)
    free(handler->handler_data);

  /*
    free this structure
   */
  if (handler!=NULL)
    free(handler);

  fprintf(stderr,"the handler at %p is freed\n",handler);

}


/*
  installs handlers in expat
 */
void
default_set_handler(object_handler* handler)
{
  XML_SetUserData(handler->parser,handler);
  XML_SetElementHandler(handler->parser,
			handler->functions.StartHandler,
			handler->functions.EndHandler);
  XML_SetCharacterDataHandler(handler->parser,
			      handler->functions.CharacterHandler);
}

/*
  pushs new handler and installs its event handlers in expat
 */
object_handler*
default_push_handler(object_handler* old_handler,
		     object_handler* new_handler)
{
  old_handler->succ=new_handler;
  new_handler->prec=old_handler;
  default_set_handler(new_handler);
  return new_handler;
}

/*
  pops this handler and installs old event handlers in expat
  returns actual event handler
 */

object_handler*
default_pop_handler(object_handler* handler)
{
  object_handler* prec;
  prec=handler->prec;
  default_set_handler(prec);
  prec->succ=NULL;
  handler->prec=NULL;
  return prec;
}

/*
  create new child
 */
void
default_StartElement_handler(void *userData,
			     const XML_Char *name,
			     const XML_Char **atts)
{
  object_handler* this_handler;
  object_handler* new_handler;

  this_handler=(object_handler*)userData;
  /* print information */
  fprintf(stderr,"New handler for %s\n",name);

  /* create new handler */
  new_handler=default_create_handler_object(this_handler->parser,name,atts);

  /* push new handler on stack */
  default_push_handler(this_handler,new_handler);
}

/*
  recieve character data from parser
 */
void
default_CharacterData_handler(void *userData,
			      const XML_Char *s,
			      int len)
{
  fprintf(stderr,"character data found\n");
}

/*
  end of a section, evaluate and test obtained data and send it to parent
 */
void
default_EndElement_handler(void *userData,
			   const XML_Char *name)
{
  object_handler* this_handler;
  this_handler=(object_handler*)userData;

  /* print information */
  fprintf(stderr,"End of handler %s\n",name);

  /* evaluate handler specific data */

  /* propagate end to parent */
  if (this_handler->prec->functions.ChildDataReciever!=NULL)
    (this_handler->prec->functions.ChildDataReciever)(this_handler);
}

/*
  recieve data from child and pop it from stack
 */
void
default_ChildData_handler(object_handler* handler)
{
  fprintf(stderr,"recieving child data at %p\n",handler);
  /* look at obtained data */

  /* pop handler from stack */
  default_pop_handler(handler);

  /* delete it */
  default_delete_handler_object(handler);
}

void
document_EndElement_handler(void *userData,
			    const XML_Char *name)
{
  object_handler* this_handler;
  this_handler=(object_handler*)userData;
  fprintf(stderr,"The document End Handler is called\n");
}

void
document_ChildData_handler(object_handler* handler)
{
  fprintf(stderr,"document_handler: recieving child data at %p\n",handler);

  /* pop handler from stack */
  default_pop_handler(handler);

  /* delete it */
  default_delete_handler_object(handler);
}

/*
  just do something...
 */
int ghmm_xml_parse(const char* filename)
{
  FILE* xml_source=NULL;
  XML_Parser my_parser;
  object_handler handler_chain_root;

  /* open file */
  xml_source=fopen(filename,"r");
  if (xml_source==NULL)
    {
      fprintf(stderr,"could not open file\n");
      return 0;
    }

  /* initialise handler chain */
  my_parser=XML_ParserCreate((const XML_Char*)NULL);
  
  /* write first structure */
  handler_chain_root.prec=NULL;
  handler_chain_root.succ=NULL;
  handler_chain_root.type=document_handler;
  handler_chain_root.parser=my_parser;
  handler_chain_root.functions.StartHandler=&default_StartElement_handler;
  handler_chain_root.functions.EndHandler=&document_EndElement_handler; /* should not be called! */
  handler_chain_root.functions.CharacterHandler=&default_CharacterData_handler;
  handler_chain_root.functions.ChildDataReciever=document_ChildData_handler;

  default_set_handler(&handler_chain_root);

  while (1)
  {
    int length;
    const int buffer_length=1000;
    void* xml_buffer;
    xml_buffer=XML_GetBuffer(my_parser,buffer_length);

    length=fread(xml_buffer,1,buffer_length,xml_source);
    if (ferror(xml_source)!=0)
      {
	fprintf(stderr,"An error occured while reading...\n");
	break;
      }

    if (XML_ParseBuffer(my_parser,length,length==0)==0)
      {
	enum XML_Error error;
	error=XML_GetErrorCode(my_parser);
	fprintf(stderr,"An error occured while parsing...\n%s\n",
		XML_ErrorString(error));	
	break;
      }

    if (length==0)
      {
	fprintf(stderr,"Ready!\n");
	break;
      }
  }

  fprintf(stderr,"deleting parser and closing file\n");
  /* delete parser */  
  XML_ParserFree(my_parser);

  /* close file */
  (void)fclose(xml_source);

  return 0;
}

#endif /* defined(HAVE_EXPAT_H) && defined(HAVE_LIBEXPAT)  */
#endif /* __experimental__ */

