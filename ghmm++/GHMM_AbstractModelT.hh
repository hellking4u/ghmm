#ifndef GHMM_ABSTRACTMODELT_H
#define GHMM_ABSTRACTMODELT_H

#include <stdio.h>
#include <vector>
#include <map>
#include <xmlio/XMLIO_Element.h>
#include <ghmm++/GHMM_Types.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Alphabet;

/** */
template<class StateType, class TransitionType> class GHMM_AbstractModelT: public XMLIO_Element {

 public:

  /** Constructor. */
  GHMM_AbstractModelT() {}
  /** Destructor. */
  ~GHMM_AbstractModelT() {;}

  /** Returns name of class. */
  const char* toString() { return "Error"; }

  /** Clean model. */
  void clean()
  {
	unsigned int i;

	for (i = 0; i < states.size(); ++i)
	  SAFE_DELETE(states[i]);
	states.clear();

	for (i = 0; i < transitions.size(); ++i)
	  SAFE_DELETE(transitions[i]);
	transitions.clear();
  }

  /**
     Tests if all standardization requirements of model are fulfilled. 
     @return 0 for succes, -1 on error
  */
  int check() { return -1; }
  
  /** Returns alphabet of model or NULL, if no such alphabet exists. */
  GHMM_Alphabet* getAlphabet() { return NULL; }

  StateType* getState(const string& id) const {
    map<string,int>::const_iterator iter;
    iter = state_by_id.find(id);
    if (iter == state_by_id.end()) {
      fprintf(stderr,"GHMM_AbstractModelTemplate::getState(id):\n");
      fprintf(stderr,"State '%s' not found in model.\n",id.c_str());
      return NULL;
    }
    
    return states[iter->second];
  }

  StateType* getState(int index) const {
    if (index >= (int) states.size()) {
      fprintf(stderr,"GHMM_AbstractModelTemplate::getState(int):\n");
      fprintf(stderr,"State no. %d does not exist. Model has %d states.\n",index,(int) states.size());
      exit(1);
    }
    return states[index];
  }

  /** Returns index of state with given id.
      @param id     id of state.
      @return index of state or -1 if no such state exists. .*/
  int getStateIndex(const string& id) const {
	  map<string,int>::const_iterator iter;
	  iter = state_by_id.find(id);
	  if (iter == state_by_id.end()) {
		return -1;
	  }
	 return iter->second;
  }
  
  /**
     Writes the model in matrix format.
     @param file: output file
  */
  void print(FILE *file) {}
  /** */
  virtual void setTransition(TransitionType* transition) {;}
  /** */
  virtual void setTransition(int start, int end, double prob) {;}
  /** */
  virtual void setTransition(const string& start, const string& end, double prob) {;}
  /** */
  void stateIDChanged(const string& old_id, const string& new_id) {
  int index = state_by_id[old_id];

  state_by_id.erase(old_id);
  state_by_id[new_id] = index;
  }

  void addStateID(const string& new_id, int index) {
	state_by_id[new_id] = index;
  }

  /** */
  vector<TransitionType*> transitions;


 protected:
  string node_tag; /* XML element-tag for a state*/
  string edge_tag; /* XML element-tag for a transition */

  void  setNodeTag( const string &tag )       { node_tag = tag; }
  void  setTransitionTag( const string &tag ) { edge_tag = tag; }

  /** Clean up c++ data structure for transitions. */
  void cleanTransitions() {
  unsigned int i;
  for (i = 0; i < transitions.size(); ++i)
    SAFE_DELETE(transitions[i]);
  transitions.clear();
  }


  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) =0;
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  virtual const int XMLIO_writeContent(XMLIO_Document& doc) = 0;

  /** */
  vector<StateType*> states;

  /** */
  map<string,int> state_by_id;
}; // class AbstractModelT<state, transition>


#ifdef HAVE_NAMESPACES
}
#endif

#endif // #ifndef 
