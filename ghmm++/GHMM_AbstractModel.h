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

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_AbstractModel;
class GHMM_State;
class GHMM_Transition;
class GHMM_Alphabet;


/** */
class GHMM_AbstractModel: public XMLIO_Element {

 public:

  /** Constructor. */
  GHMM_AbstractModel();
  /** Destructor. */
  virtual ~GHMM_AbstractModel();

  /** Returns name of class. */
  virtual const char* toString() const;

  /** */
  void addStateID(const string& id, int index);
  /**
     Tests if all standardization requirements of model are fulfilled. 
     @return 0 for succes, -1 on error
  */
  virtual int check() const;
  /** Clean model. */
  virtual void clean();
  /** Returns alphabet of model or NULL, if no such alphabet exists. */
  virtual GHMM_Alphabet* getAlphabet() const;
  /** Returns model type. */
  virtual GHMM_ModelType getModelType() const;
  /** */
  virtual int getNumberOfTransitionMatrices() const;
  /* Returns state with given index. */
  GHMM_State* getState(int index) const;
  /* Returns state with given id. */
  GHMM_State* getState(const string& id) const;
  /** Returns index of state with given id.
      @param id     id of state.
      @return index of state or -1 if no such state exists. .*/
  int getStateIndex(const string& id) const;
  /**
     Writes the model in matrix format.
     @param file: output file
  */
  virtual void print(FILE *file) const;
  /** */
  void setTransition(GHMM_Transition* transition);
  /** */
  void setTransition(int start, int end, double prob);
  /** */
  void setTransition(const string& start, const string& end, double prob);
  /** */
  void stateIDChanged(const string& old_id, const string& new_id);

  /** */
  vector<GHMM_Transition*> transitions;


 protected:
  string node_tag; /* XML element-tag for a state*/
  string edge_tag; /* XML element-tag for a transition */

  void  setNodeTag( const string &tag );
  void  setTransitionTag( const string &tag );

  /** Clean up c++ data structure for transitions. */
  void cleanTransitions();
  /** Build up c++ data structure for transitions. */
  void createTransitions();
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  virtual const int XMLIO_writeContent(XMLIO_Document& doc);

  /** */
  vector<GHMM_State*> states;
  /** */
  map<string,int> state_by_id;
};



#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_ABSTRACTMODEL_H */
