/*
 * created: 05 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#ifndef _GHMM_ABSTRACTMODEL_H
#define _GHMM_ABSTRACTMODEL_H 1

#include <stdio.h>
#include <vector>
#include <xmlio/XMLIO_Element.h>
#include <ghmm++/GHMM_Types.h>

#include <ghmm++/GHMM_AbstractModelT.hh>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_State;
class GHMM_Transition;

class GHMM_AbstractModel : public GHMM_AbstractModelT<GHMM_State, GHMM_Transition>
{

public:
  GHMM_AbstractModel();

	/** Destructor. */
  ~GHMM_AbstractModel();
  
  /** Returns name of class. */
  const char* toString() const;
  
  /**
     Tests if all standardization requirements of model are fulfilled. 
     @return 0 for succes, -1 on error
  */
  int check() const;
  
  
  /** Returns alphabet of model or NULL, if no such alphabet exists. */
  GHMM_Alphabet* getAlphabet() const;
  /** Returns model type. */
  virtual GHMM_ModelType getModelType() const = 0;
  /** */
  virtual int getNumberOfTransitionMatrices() const = 0;

  /** */
  virtual void *get_cmodel() const = 0;

  /* Returns state with given index. */
  //GHMM_State* getState(int index) const;
  /* Returns state with given id. */
  //GHMM_State* getState(const string& id) const;

   /**
     Writes the model in matrix format.
     @param file: output file
  */
  void print(FILE *file) const;
  /** */
  void setTransition(GHMM_Transition* transition);
  /** */
  void setTransition(int start, int end, double prob);
  /** */
  void setTransition(const string& start, const string& end, double prob);
  
protected:
  
  void createTransitions();

	  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  const int XMLIO_writeContent(XMLIO_Document& doc);

  void setNodeTag( const string &tag );

  void setTransitionTag( const string &tag );

}; // class AbstractModel


#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_ABSTRACTMODEL_H */
