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
#include <xmlio/XMLIO_AttributedElement.h>
#include <xmlio/XMLIO_ArrayElement.h>


#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
     collect everything with weight
   */
  class Emission: public XMLIO_ContentElementArrayElement<double,XMLIO_AttributedElement>
    {
    public:
      /** Poor Constructor. */
      Emission();

      /** with attributes State and fix */
      Emission(const string& tag, XMLIO_Attributes &attributes);
      /** */
      virtual const char* toString() const;
      /** name of the state */
      string state;
      /** if fix==1 the emissions of this state is not trained */
      int fix;
    };

  /**
   */
  class Emissions: public XMLIO_ElementArrayElement<Emission>
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
  
