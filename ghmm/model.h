/*******************************************************************************
  author       : Bernd Wichern
  filename     : ghmm/ghmm/model.h
  created      : TIME: 10:47:27     DATE: Fri 19. December 1997
  $Id$

__copyright__

*******************************************************************************/


#ifndef MODEL_H
#define MODEL_H

/**@name HMM-Modell */
/*@{ (Doc++-Group: model) */

/** @name state
    Grundlegende Struktur, in der alle zu einem Zustand gehoerigen Parameter
    gespeichert werden.
*/
struct state{
  /** Anfangswahrscheinlichkeit */ 
  double pi;  
  /** Ausgabewahrscheinlichkeit */
  double *b;
  /** ID's der Nachfolge States */ 
  int *out_id;  
  /** ID's der Vorgaenger States */    
  int *in_id;
  /** Nachfolge Uebergangswahrscheinlichkeiten */       
  double *out_a; 
  /** Vorgaenger Uebergangswahrscheinlichkeiten */   
  double *in_a;
  /** Anzahl Nachfolge States */     
  int out_states; 
  /** Anzahl Vorgaenger States */
  int in_states;  
  /** falls fix == 1 --> b bleibt fix, wird also nicht trainiert */
  int fix;
};
typedef struct state state;

/** @name struct model
    Gesamt HMM-Modell. Enth\"alt also alle Parameter, die ein HMM festlegen.
*/
struct model{
  /** Zahl der Zustaende */
  int N;
  /** Zahl der Ausgabewerte */   
  int M;   
  /** Vektor der Zustaende*/
  state *s; 
  /** Prior fuer die a priori Wahrscheinlichkeit des Modells.
   Kein prior definiert: Wert auf -1 gesetzt */
  double prior;
};
typedef struct model model;


/** @name model\_direct
    Modellparameter in Matrixform fuer hmm\_input
*/
struct model_direct{
  /** Zahl der Zustaende */
  int N;
  /** Zahl der Ausgabewerte */  
  int M;
  /** Prior fuer die a priori Wahrscheinlichkeit des Modells.
   Kein prior definiert: Wert auf -1 gesetzt */
  double prior;
  /** Uebergangsmatrix */
  double **A;
  /** Ausgabematrix */
  double **B;
  /** Anfangsmatrix */
  double * Pi;
  /** Vektor zur Kennzeichnung von States, deren Ausgabeparameter nicht 
      trainiert werden sollen. Default ist 0 fuer alle states */
  int *fix_state;
};
typedef struct model_direct model_direct;

/** @name hmm_check\_t
  */
struct hmm_check_t{
  /** */
  int r_a;
  /** */
  int c_a;
  /** */
  int r_b;
  /** */
  int c_b;
  /** */
  int len_pi;
  int len_fix;
};
typedef struct hmm_check_t hmm_check_t;

/*
  Wichtig: Include von sequence.h zur Vermeidung von Fehlern beim Uebersetzen
  erst an dieser Stelle. (sequence.h und model.h includen sich gegenseitig)
*/
#include "sequence.h" 
#include "scanner.h"

/** Speicher eines Modells freir\"aumen..
    @return Erfolgsstatus: 0 succes; -1 error
    @param mo:  ptr. auf ptr auf Modell */
int     model_free(model **mo);

/**
   Einlesen einer Ascii-Datei zur Initialisierung der Modelle. Alloc vom 
   Speicher f\"ur die Modelle erfolgt hier.
   @return Vektor von Pointern der eingelesenen Modelle
   @param filename:   Name der einzulesenden Ascii-Datei
   @param mo_number:  Zahl der eingelesenen Modelle */
model** model_read(char *filename, int *mo_number);
/**
   Einlesen eines Modells bei expliziter Angabe aller relevanten 
   Modellparameter in Matrixform. Alloc vom Speicher f\"ur die Modelle erfolgt hier.
   @return Pointer auf eingelesenes Modell
   @param s:       scanner
   @param multip:  eingelesene Multiplizitaet; gibt Anzahl der Kopien an, die vom 
   eingelesenen Modell erzeugt werden sollen */
model*  model_direct_read(scanner_t *s, int *multip);

/**
   Erzeugung von einfachen left-right Modellen nach vorgegebenen Sequenzen.
   Daf\"ur wird f\"ur jede Sequenz "model_generate_from_sequence" aufgerufen.
   Die Sequenzen werden direkt aus dem Ascii File eingelesen und nach 
   Verlassen der Funktion wieder verworfen.
   @return Vektor von Modellen
   @param s:          scanner
   @param new_models: Anzahl der erzeugten Modelle */

model **model_from_sequence_ascii(scanner_t *s, long *mo_number);


/** 
    Aufgabe siehe: model_from_sequence_ascii. Unterschied: die Sequenzen werden
    nicht aus Ascii File eingelesen, sondern liegen bereits als Struktur vor.
    @return Vektor von Modellen
    @param s:          scanner
    @param new_models: Anzahl der erzeugten Modelle */

model **model_from_sequence(sequence_t *sq, long *mo_number);

/**
   Kopieren eines vorgegebenen Modells. Speicher alloc geschieht in model\_copy.
   @return Kopie des Modells
   @param mo:  zu kopierendes Modell */
model*  model_copy(const model *mo);
/**
   Pr\"ufen, ob alle Normierungsbedingungen vom Modell erf\"ullt werden. (Summen 
   der versch. Wahrscheinlichkeiten = 1?)
   @return Erfolgsstatus: 0 succes; -1 error
   @param mo:  zu ueberpruefendes Modell */
int     model_check(const model* mo);
/**
   Check Anzahl States und Anzahl Ausgabewerte in den Modellen auf
   Uebereinstimmung..
   @return Erfolgsstatus: 0 succes; -1 error
   @param mo:           Vektor der Modelle
   @param model_number: Anzahl der Modelle */
int     model_check_compatibility(model **mo, int model_number);
/**
   Erzeugt ein Modell, das die uebergebene Sequenz mit P = 1 ausgibt. 
   Striktes Left-Right Model, pro Element der
   Sequenz ein Zustand, in dem mit P = 1 das gegebene Symbol ausgegeben wird.
   Das Modell erhaelt speziellen "final State", der keine out-states hat.
   @return         Pointer auf das erzeugte Modell (Allozierung hier)
   @param seq      Vorgabesequenz
   @param seq_len  Laenge der Vorgabesequenz
   @param anz_symb Zahl der Symbole der Vorgabesequenz
*/
model*  model_generate_from_sequence(const int *seq, int seq_len, 
				     int anz_symb);

/** 
    Erzeugt zu einem gegebenen Modell zuf\"allige Sequenzen.
    Speicher f\"ur die Sequenzen und den L\"angenvektor wird in der Fkt. selbst
    bereitgestellt.
    Die L\"ange der Sequenzen kann global vorgegeben werden (global_len > 0)
    oder in der Funktion dadurch bestimmt werden, dass ein "final state"
    erreicht wird (ein state, dessen Ausgangswahrscheinlichkeiten =0 sind).
    Besitzt das Modell keinen Final State, so werden Sequenzen der L\"ange 
    MAX_SEQ_LEN generiert. 
    @return            Pointer auf ein Feld von Sequenzen (Allozierung)
    @param mo          vorgegebenes Modell
    @param seed        Initialisierungsvariable des Zufallsgenerators (int).
                       Wird 0 übergeben, wird der RNG nicht ge-seeded 
    @param global_len  gewuenschte Sequenzlaenge (=0: autom. ueber final state)
    @param seq_number  gewuenschte Anzahl von Sequenzen
*/
sequence_t *model_generate_sequences(model* mo, int seed, int global_len,
				     long seq_number);


/**
   Berechnung der Summe der log( P ( O | lambda ) ).
   Sequenzen, die vom Modell nicht dargestellt werden koennen, werden 
   vernachlaessigt.
   @return    log(P)
   @param mo vorgegebenes Modell
   @param sq vorgegebene Sequenzen       
*/
double model_likelihood(model *mo, sequence_t *sq);
/**
   Schreiben eines Modells in Matrixform..
   @param file: Ausgabedatei
   @param mo:   zu schreibendes Modell
*/
void model_print(FILE *file, model *mo); 
/**
   Schreiben der Uebergangsmatrix eines Modells..
   @param file: Ausgabedatei
   @param mo:   zu schreibendes Modell
   @param tab:  zur Formatierung: fuehrende Tabs
   @param separator: Formatierung: Trennzeichen der Spalten
   @param ending:    Formatierung: Endzeichen der Zeile  
*/
void model_A_print(FILE *file, model *mo, char *tab, char *separator, 
		   char *ending);
/**
   Schreiben der Ausgabematrix eines Modells..
   @param file: Ausgabedatei
   @param mo:   zu schreibendes Modell
   @param tab:  zur Formatierung: fuehrende Tabs
   @param separator: Formatierung: Trennzeichen der Spalten
   @param ending:    Formatierung: Endzeichen der Zeile  
*/
void model_B_print(FILE *file, model *mo, char *tab, char *separator, 
		   char *ending);
/**
   Schreiben des Anfangsbelegungsvektors eines Modells..
   @param file: Ausgabedatei
   @param mo:   zu schreibendes Modell
   @param tab:  zur Formatierung: fuehrende Tabs
   @param separator: Formatierung: Trennzeichen der Spalten
   @param ending:    Formatierung: Endzeichen der Zeile  
*/
void model_Pi_print(FILE *file, model *mo, char *tab, char *separator, 
		    char *ending);

void model_fix_print(FILE *file, model *mo, char *tab, char *separator, 
		     char *ending);
/**
   Schreiben der transponierten Uebergangsmatrix eines Modells..
   @param file: Ausgabedatei
   @param mo:   zu schreibendes Modell
   @param tab:  zur Formatierung: fuehrende Tabs
   @param separator: Formatierung: Trennzeichen der Spalten
   @param ending:    Formatierung: Endzeichen der Zeile  
*/
void model_A_print_transp(FILE *file, model *mo, char *tab, char *separator, 
			  char *ending);
/**
   Schreiben der transponierten Ausgabematrix eines Modells..
   @param file: Ausgabedatei
   @param mo:   zu schreibendes Modell
   @param tab:  zur Formatierung: fuehrende Tabs
   @param separator: Formatierung: Trennzeichen der Spalten
   @param ending:    Formatierung: Endzeichen der Zeile  
*/
void model_B_print_transp(FILE *file, model *mo, char *tab, char *separator, 
			  char *ending);
/**
   Schreiben des transponierten Anfangsbelegungsvektors eines Modells..
   @param file: Ausgabedatei
   @param mo:   zu schreibendes Modell
   @param tab:  zur Formatierung: fuehrende Tabs
   @param separator: Formatierung: Trennzeichen der Spalten
   @param ending:    Formatierung: Endzeichen der Zeile  
   */
void model_Pi_print_transp(FILE *file, model *mo, char *tab, char *ending);
/**
   Schreiben eines HMM-Modells in Matrixformat. 
   Input: model\_direct (heisst: Parameter als Matrix gespeichert.)
   @param file:   Ausgabedatei
   @param mo_d:   zu schreibendes Modell in Matrixformat
   @param multip: Modell wird in multip Kopien geschrieben
*/
void model_direct_print(FILE *file, model_direct *mo_d, int multip);
/** 
    Etwas unuebersichtliche Parameterausgabe sortiert nach den States 
    des Modells..
    @param file: Ausgabedatei
    @param mo:   zu schreibendes Modell
*/
void model_states_print(FILE *file, model *mo); 

/** Alle Pointer des Modell freigeben und auf NULL setzen, Variablen 
    auf Null setzen..
    @param mo_d  HMM-Modell-Struktur (\Ref{struct model_direct})
    @param check Check-Struktur (\Ref{struct hmm_check_t})
*/
void model_direct_clean(model_direct *mo_d, hmm_check_t *check); 

/** Teste Kompatibiliatet der Modell-Komponenten..
    @return ??
    @param mo_d  HMM-Modell-Struktur  (\Ref{struct model_direct})
    @param check Check-Struktur  (\Ref{struct hmm_check_t})
*/
int model_direct_check_data(model_direct *mo_d, hmm_check_t *check); 

/** Computes probabilistic distance of two models
    @return the distance
    @param m0  model used for generating random output
    @param m  model to compare with
    @param maxT  maximum output length (for HMMs with absorbing states multiple
                 sequences with a toal langth of at least maxT will be 
		 generated)
    @param symmetric  flag, whether to symmetrize distance (not implemented yet)
    @param verbose  flag, whether to monitor distance in 40 steps. 
                    Prints to stdout (yuk!)
*/
double model_prob_distance(model *m0, model *m, int maxT, int symmetric, int verbose);






#endif

/*@} (Doc++-Group: model) */










