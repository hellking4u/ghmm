/*-----------------------------------------------------------------------------
  author       : Alexander Schliep
  filename     : ghmm.h
  created      : TIME: 21:44:45     DATE: Sun 14. November 1999
  $Id$

__copyright__

------------------------------------------------------------------------------*/

#ifndef GHMM_H
#define GHMM_H

#ifdef __cplusplus
extern "C" {
#endif
/**@name GHMM-Globals */
/*@{ (Doc++-Group: globals) */
#ifdef __cplusplus
}
#endif



/** @name type_constants
    Constants giving model variations */
/** Model is a left-right */
#define kNotSpecified (0)
#define kLeftRight (1)
/** Model contains silent states (i.e., states without emissions) */
#define kSilentStates (1 << 2)
/** Model has states with tied emission probabilities */
#define kTiedEmissions (1 << 3)
#define kUntied -1 

/** Model has states emission probabilities conditioned on previous orders */
#define kHigherOrderEmissions (1 << 4);

#define kHasBackgroundDistributions (1 << 5);

#define kNoBackgroundDistribution -1


#endif

/*@} (Doc++-Group: GHMM-Globals) */

