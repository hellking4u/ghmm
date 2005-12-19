
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


xmlChar *ConvertInput(const char *in, const char *encoding);


int writeIntAttribute(xmlTextWriterPtr writer, const char * name, int value) {
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
void writeHMMDocument(fileData_s * f, const char *file) {
#define CUR_PROC "writeHMMDocument"
  int rc;
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
    printf("testXmlwriterDoc: Error at xmlTextWriterStartElement\n");
    return;
  }

  if (writeIntAttribute(writer, "noComponents", f->noModels)) {
    GHMM_LOG_QUEUED(LERROR);
    goto STOP;
  }


  /* Here we could close the elements ORDER and EXAMPLE using the
   * function xmlTextWriterEndElement, but since we do not want to
   * write any other elements, we simply call xmlTextWriterEndDocument,
   * which will do all the work. */
  rc = xmlTextWriterEndDocument(writer);
  if (rc < 0) {
    GHMM_LOG(LERROR, "Error at xmlTextWriterEndDocument");
    goto STOP;
  }

  xmlSaveFileEnc(file, doc, MY_ENCODING);

STOP:
  xmlFreeTextWriter(writer);
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
xmlChar *
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
