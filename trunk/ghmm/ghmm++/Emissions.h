/*
  author: Achim Gaedke
  created: 13. Juli 2001
  file: xmlio/examples/ghmm++/Emissions.h
  $Id$
 */

#ifndef GHMMPP_EMISSIONS_H
#define GHMMPP_EMISSIONS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <xmlio/XMLIO_AttributedObject.h>
#include <xmlio/XMLIO_ArrayObject.h>


#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
     collect everything with weight
   */
  class Emission: public XMLIO_ContentElementArrayObject<double,XMLIO_AttributedObject>
    {
    public:
      /** Poor Constructor. */
      Emission();

      /** with attributes State and fix */
      Emission(const string& tag, XMLIO_Attributes &attributes);
      /** */
      virtual const char* toString() const;

    protected:
      /** */
      string state;
      /** */
      int fix;
    };

  /**
   */
  class Emissions: public XMLIO_ElementArrayObject<Emission>
    {
    public:
      /** Poor constructor. */
      Emissions();

      /** sets the name of collected Elements to Emission */
      Emissions(const string& tag, XMLIO_Attributes &attributes);

      /** */
      virtual const char* toString() const;

    };

#ifdef HAVE_NAMESPACES
}
#endif
  
#endif /* GHMMPP_EMISSIONS_H */
  
