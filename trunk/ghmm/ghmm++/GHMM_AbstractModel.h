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

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_AbstractModel;
class GHMM_State;
class GHMM_Transition;

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
  void setTransition(GHMM_Transition* transition);
  /** */
  void setTransition(int start, int end, double prob);
  /**
     Tests if all standardization requirements of model are fulfilled. 
     @return 0 for succes, -1 on error
  */
  virtual int check() const;
  /** Clean model. */
  virtual void clean();
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
  vector<GHMM_Transition*> transitions;


 protected:

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
