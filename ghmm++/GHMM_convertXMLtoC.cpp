/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#include <stdio.h>

#include "ghmm++/GHMM.h"
#include "ghmm++/GHMM_ContinuousModel.h"
#include "ghmm++/GHMM_GMLDoc.h"
#include "ghmm++/GHMM_convertXMLtoC.h"
// #include "ghmm++/templatetest.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


//myvector<double> Vecdouble;

/** This function is compiled as a C-function ( no C++ "name-mangled" )
    so that we can read data from an XML-file into a C "struct",
    ( struct model and struct smodel )
    @return Return a pointer to struct model_t which provides
    access to a C struct whose type must be forced at run-time by a user.
    Its type-tag is stored in a member variable "model_id" of model_t.
*/
model_t *graphmldoc_cwrapper(char *filename)
{
  /*   GHMM_GraphMLDoc    doc; */
  GHMM_GraphMLDoc    doc;
  // GHMM_GMLDiscreteModel *dmo;
  GHMM_SWDiscreteModel *dmo;
  GHMM_ContinuousModel *cmo;

  model_t *model_pt;
  void    *tmpvoid;

  doc.open(filename, "r");
  doc.readDocument();
  doc.close();

  if ( (tmpvoid = calloc( 1, sizeof(model_t))) != NULL )
    {
      model_pt = (model_t*) tmpvoid;
    } 
  else
    {
      fprintf(stderr, "Cannot allocate memory\n");
      model_pt = NULL; /* cannot allocate memory*/
      return model_pt;
    }

  if ( doc.getDiscreteModel() != NULL )
    {
      model_pt->model_id = DISCRETE;
      dmo = doc.getDiscreteModel();
      //
      // Make a new copy so that we can free the pointer to a C++ object. 
      //
      model_pt->model_pt = (void*) sdmodel_copy((const sdmodel*) dmo->c_model);
      delete dmo;
    }
  else
    {
      model_pt->model_id = CONTINUOUS;
      cmo = doc.getContinuousModel();
      //
      // Make a new copy so that we can free the pointer to a C++ object. 
      //
      model_pt->model_pt = (void*) smodel_copy((const smodel*) cmo->c_model);
      delete dmo;
    }

  return model_pt;
}
