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
default_create_handler_object(document_handler* dh,
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
  new_handler->type=default_object;
  new_handler->dh=dh;
  new_handler->handler_data=NULL;
  new_handler->functions.StartHandler=&default_StartElement_handler;
  new_handler->functions.EndHandler=&default_EndElement_handler;
  new_handler->functions.CharacterHandler=&default_CharacterData_handler;
  new_handler->functions.ChildDataReciever=&default_ChildData_handler;

  /*
    evaluate type, name and attributes
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
  create new child
 */
void
default_StartElement_handler(void *userData,
			     const XML_Char *name,
			     const XML_Char **atts)
{
  document_handler* my_dh;
  object_handler* new_handler;
  my_dh=((object_handler*)userData)->dh;

  /* print information */
  fprintf(stderr,"New handler for %s\n",name);

  /* create new handler */
  new_handler=create_object_handler(my_dh,
				    default_object,
				    name,
				    atts);

  /* push new handler on stack */
  document_handler_push_object_handler(my_dh,new_handler);
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
  document_handler_pop_object_handler(handler->dh);

  /* delete it */
  delete_object_handler(handler);
}

/***************************************************************/
/* abstract object_handler functions */
object_handler*
create_object_handler(document_handler* dh,
		      const object_handler_type type,
		      const XML_Char *name,
		      const XML_Char **atts)
{
  /* nothing other implemented ! */
  return default_create_handler_object(dh,name,atts);
}

void
delete_object_handler(object_handler* handler)
{
  default_delete_handler_object(handler);
}


/***************************************************************/

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
  (void)document_handler_pop_object_handler(handler->dh);

  /* delete it */
  delete_object_handler(handler);
}


/***************************************************************/


/* create document handler */

document_handler*
create_document_handler(const char* filename)
{
  document_handler* my_dh;
  my_dh=(document_handler*)malloc(sizeof(document_handler));
  if (my_dh==NULL)
    {
      fprintf(stderr,"Can not allocate document_handler\n");
      return NULL;
    }

  /* open file */
  my_dh->xml_file=fopen(filename,"r");
  if (my_dh->xml_file==NULL)
    {
      fprintf(stderr,"could not open file\n");
      free(my_dh);
      return NULL;
    }

  /* initialise handler chain */
  my_dh->parser=XML_ParserCreate((const XML_Char*)NULL);
  if (my_dh->parser==NULL)
    {
      fprintf(stderr,"could not create parser!\n");
      fclose(my_dh->xml_file);
      free(my_dh);
      return NULL;
    }

  /* write first structure */
  my_dh->root.prec=NULL;
  my_dh->root.succ=NULL;
  my_dh->last=&(my_dh->root);
  my_dh->root.type=document_object;
  my_dh->root.functions.StartHandler=&default_StartElement_handler;
  my_dh->root.functions.EndHandler=&document_EndElement_handler; /* should not be called! */
  my_dh->root.functions.CharacterHandler=&default_CharacterData_handler;
  my_dh->root.functions.ChildDataReciever=document_ChildData_handler;

  document_handler_set_handler(&(my_dh->root));

  return my_dh;
}

/*
  installs handlers in expat
 */
void
document_handler_set_handler(object_handler* handler)
{
  XML_Parser* parser;
  parser=handler->dh->parser;
  XML_SetUserData(parser,handler);
  XML_SetElementHandler(parser,
			handler->functions.StartHandler,
			handler->functions.EndHandler);
  XML_SetCharacterDataHandler(parser,
			      handler->functions.CharacterHandler);
}

/*
  pushs new handler on stack and installs its event handlers in expat
 */
object_handler*
document_handler_push_object_handler(document_handler* dh,
		     object_handler* new_handler)
{
  dh->last->succ=new_handler;
  new_handler->prec=dh->last;
  new_handler->succ=NULL;
  dh->last=new_handler;
  document_handler_set_handler(new_handler);
  return new_handler;
}

/*
  pops this handler and installs old event handlers in expat
  returns removed handler
 */
object_handler*
document_handler_pop_object_handler(document_handler* dh)
{
  object_handler* last;
  last=dh->last;
  dh->last=last->prec;
  dh->last->prec=NULL;
  last->prec=last->succ=NULL;
  document_handler_set_handler(dh->last);
  return last;
}

int
parse_document(document_handler* dh)
{
  const int buffer_length=1000;
  while (1)
    {
      int length;
      void* xml_buffer;
      xml_buffer=XML_GetBuffer(dh->parser,buffer_length);
      
      length=fread(xml_buffer,1,buffer_length,dh->xml_file);
      if (ferror(dh->xml_file)!=0)
	{
	  fprintf(stderr,"An error occured while reading...\n");
	  return 0;
	}
      
      if (XML_ParseBuffer(dh->parser,length,length==0)==0)
	{
	  enum XML_Error error;
	  error=XML_GetErrorCode(dh->parser);
	  fprintf(stderr,"An error occured while parsing...\n%s\n",
		  XML_ErrorString(error));	
	  return 0;
	}
      
      if (length==0)
	{
	  fprintf(stderr,"Ready!\n");
	  break;
	}
    }
  return 1;
}

/* destroy document handler */

void
delete_document_handler(document_handler* dh)
{
  /* rebuild stack */
  /* delete all, except root*/
  while (dh->last->prec!=NULL)
    {
      object_handler* old_object_handler;
      old_object_handler=document_handler_pop_object_handler(dh);
      delete_object_handler(old_object_handler);
    }

  XML_ParserFree(dh->xml_file);
  fclose(dh->xml_file);
  free(dh);
}

#endif /* defined(HAVE_EXPAT_H) && defined(HAVE_LIBEXPAT)  */
#endif /* __experimental__ */


