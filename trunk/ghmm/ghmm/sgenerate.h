/*------------------------------------------------------------------------------
  author       : Bernhard Knab
  filename     : ghmm/ghmm/sgenerate.h
  created      : TIME: 09:33:23     DATE: Tue 16. November 1999
  $Id$

__copyright__

------------------------------------------------------------------------------*/
#ifndef SGENERATE_H
#define SGENERATE_H


#include <ghmm/smodel.h>
#include <ghmm/sequence.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
   @name generation and extention of sequences from shmm
*/

/*@{
 */

/**
 Kopiert von cgenerate.h
*/

typedef enum {
  viterbi_viterbi,
  viterbi_all,
  all_viterbi,
  all_all
} sgeneration_mode_t;

/**
   Verlaengert Teilsequenzen bei gegebenem Modell. Verschiedene Moeglichkeiten
   sind vorgesehen (nur Viterbi oder alle Pfade beruecksichtigen kombiniert
   mit Sequenzanfang und Sequenzende)
   Steuerung ueber mode:
   0 = viterbi\_viterbi, 
   1 = viterbi\_all, 
   2 = all\_viterbi, 
   3 = all\_all
   (zunaechst nur all\_all moeglich)
   Die generierten Sequenzen werden zurueckgegeben.
   @return Pointer auf ein Feld von Gesamtsequenzen (vorgegebene 
   Anfangssequenz und generierte Endsequenz)
   @param smo         vorgegebenes Modell
   @param sqd_short   Feld von Anfangssequenzen
   @param seed        Initialisierungsvariable des Zufallsgenerators (int)
   @param global_len  gewuenschte Sequenzlaenge (=0: autom. ueber final state)
   @param mode        Steuerung der Generierung
 */
sequence_d_t *sgenerate_extensions(smodel *smo, sequence_d_t *sqd_short, 
				   int seed, int global_len,
				   sgeneration_mode_t mode);


/** Verlaengern einer einzelnen Anfangsequenz. Sonst gleiche Funktionalitaet
    wie sgenerate_extensions
*/
double *sgenerate_single_ext(smodel *smo, double *O, const int len, 
			     int *new_len, double **alpha,
			     sgeneration_mode_t mode);


/** generate a single next value bases on a trained model and on a seq und
   to length "len"
*/
double sgenerate_next_value(smodel *smo, double *O, const int len);

/*@} sgenerate section */

#ifdef __cplusplus
}
#endif


#endif /* SGENERATE_H */
