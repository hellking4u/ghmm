/*-----------------------------------------------------------------------------
  author       : Bernhard Knab
  filename     : /homes/hmm/wichern/hmm/src/smodel.h
  created      : TIME: 21:44:45     DATE: Sun 14. November 1999
  last-modified: TIME: 14:19:17     DATE: Wed 07. February 2001
------------------------------------------------------------------------------*/

#ifndef SMODEL_H
#define SMODEL_H

#include "const.h"
#include "scanner.h"

/**@name SHMM-Modell */
/*@{ (Doc++-Group: smodel) */

/** @name GET\_DELTA\_OSC
    Globales define
*/
#define GET_DELTA_OSC(c, class) ((c) == (class) ? 1 : 0)

/* @name smodel\_get\_osc
    Bestimmung der osum-Klasse

    int smodel_get_osc(double osum);
*/


/** @name sstate
    Grundlegende Struktur, in der alle zu einem Zustand gehoerigen Parameter
    gespeichert werden.
*/
struct sstate{
  /** Anfangswahrscheinlichkeit */ 
  double pi;  
  /** ID's der Nachfolge States */ 
  int *out_id;  
  /** ID's der Vorgaenger States */    
  int *in_id;
  /** Nachfolge Uebergangswahrscheinlichkeiten */       
  double **out_a; 
  /** Vorgaenger Uebergangswahrscheinlichkeiten */   
  double **in_a;
  /** Anzahl Nachfolge States */     
  int out_states; 
  /** Anzahl Vorgaenger States */
  int in_states;    
  /** Gewichtsvektor. Vektor der Gewichte der versch. Verteilungen */
  double *c;
  /** Mittelwertvektor. Vektor der Mittelwerte der versch. Verteilungen */
  double *mue;
  /** Varianzvektor. 
      Vektor der Varianzen der versch. Verteilungen */
  double *u;
  /** falls fix == 1 --> c, mue, u bleiben fix, werden also nicht trainiert */
  int fix;
};
typedef struct sstate sstate;

/** @name smodel
    continous HMM
*/
struct smodel{
  /** Zahl der Zustaende */
  int N;
  /** Zahl der Verteilungen pro Zustand */   
  int M;
  /** Flag fuer die zu verwendende Dichtefunktion */
  density_t density;
  /** Prior fuer die a priori Wahrscheinlichkeit des Modells.
   Kein prior definiert: Wert auf -1 gesetzt */
  double prior;
  /** Vektor der Zustaende */
  sstate *s; 
};
typedef struct smodel smodel;

/** @name smodel\_direct
    Modellparameter in Matrixform fuer hmm\_input
*/
struct smodel_direct{
  /** Zahl der Zustaende */
  int N;
  /** Zahl der Verteilungen */  
  int M;
  /** Flag fuer die zu verwendende Dichtefunktion */
  density_t density;
  /** Prior fuer die a priori Wahrscheinlichkeit des Modells.
   Soll kein Prior verwendet werden, Wert auf -1 setzen. */
  double prior;
  /** Anfangsmatrix */
  double * Pi;
  /** 3-dim. Uebergangsmatrix (BEACHTE: O-Sum-Klasse k als erster Index) */
  double ***A;
  /** Matrix der Verteilungsgewichte */
  double **C;
  /** Matrix der Mittelwerte (NxM) */
  double **Mue;
  /** Matrix der Varianzen (NxM) */
  double **U;
  /** Vektor zur Kennzeichnung von States, deren Ausgabeparameter nicht 
      trainiert werden sollen. Default ist 0 fuer alle states */
  int *fix_state;
};
typedef struct smodel_direct smodel_direct;


/** @name shmm_check\_t
  */
struct shmm_check_t{
  /** */
  int r_a;
  /** */
  int c_a;
  /** */
  int r_c;
  /** */
  int c_c;
  /** */
  int r_mue;
  /** */
  int c_mue;
  /** */
  int r_u;
  /** */
  int c_u;
  /** */
  int len_pi;
  int len_fix;
};
typedef struct shmm_check_t shmm_check_t;

/*
  Wichtig: Include von sequence.h zur Vermeidung von Fehlern beim Uebersetzen
  erst an dieser Stelle. (sequence.h und smodel.h includen sich gegenseitig) 
*/

#include "sequence.h"  


/** Speicher eines Modells freiraeumen..
    @return ?? 
    @param smo  ptr. auf ptr auf Modell */
int     smodel_free(smodel **smo);

/**
   Einlesen einer Ascii-Datei zur Initialisierung der Modelle..
   @return Vektor der eingelesenen Modelle
   @param filename    Name der einzulesenden Ascii-Datei
   @param smo_number  Zahl der eingelesenen Modelle */
smodel** smodel_read(char *filename, int *smo_number);

/**
   Einlesen eines Modells bei expliziter Angabe aller relevanten 
   Modellparameter in Matrixform..
   @return ??
   @param s        scanner
   @param multip   eingelesene Multiplizitaet; gibt Anzahl der Kopien an, 
                   die vom eingelesenen Modell erzeugt werden sollen */

smodel*  smodel_direct_read(scanner_t *s, int *multip);

/**
   Kopieren eines vorgegebenen Modells. Speicher alloc geschieht in 
   smodel\_copy.
   @return ??
   @param smo   zu kopierendes Modell */
smodel*  smodel_copy(const smodel *smo);

/**
   Pruefen, ob alle Normierungsbedingungen vom Modell erfuellt werden. (Summen 
   der versch. Wahrscheinlichkeiten = 1?)
   @return ??
   @param smo   zu ueberpruefendes Modell */
int     smodel_check(const smodel* smo);

/**
   Check Anzahl States und Anzahl Ausgabewerte in den Modellen auf
   Uebereinstimmung..
   @return ??
   @param smo             Vektor der Modelle
   @param smodel_number   Anzahl der Modelle */
int     smodel_check_compatibility(smodel **smo, int smodel_number);

/**
   Erzeugt eine Zufallszahl der m-ten Dichtefunktion im Zustand state 
   des Modells, entsprechend der Art der Dichtefunktion des Modells.
   Greift dabei auf Fkt. in randvar zurueck.
   @return                erzeugter Zufallszahlwert
   @param smo             smodel
   @param state           Zustand
   @param m               Index der Dichtefunktion im Zustand */
double smodel_get_random_var(smodel *smo, int state, int m);

/**
   Berechnung der Summe der log( P ( O|lambda ) ).
   Sequenzen, die vom Modell nicht dargestellt werden koennen, werden 
   nicht mitsummiert.
   @return       n/-1 wobei  n = Anzahl der eingegangenen Seq./ -1 = Fehler
   @param smo    vorgegebenes Modell
   @param sq     vorgegebene Sequenzen       
   @param log\_p  gesuchte Likelihood (Rückgabe)
*/
int smodel_likelihood(smodel *smo, sequence_d_t *sqd, double *log_p);

/**
   Schreiben eines Modells in Matrixform..
   @param file   Ausgabedatei
   @param smo    zu schreibendes Modell
*/
void smodel_print(FILE *file, smodel *smo); 

/**
   Schreiben eines Modells in Matrixform MIT NUR EINER MATRIX A (= Ak_0).
   @param file   Ausgabedatei
   @param smo    zu schreibendes Modell
*/
void smodel_print_oneA(FILE *file, smodel *smo);

/**
   Schreiben der Uebergangsmatrix einer O-Sum-Klasse eines Modells..
   @param file       Ausgabedatei
   @param smo        zu schreibendes Modell
   @param k          O-Sum-Klasse k
   @param tab        zur Formatierung: fuehrende Tabs
   @param separator  Formatierung: Trennzeichen der Spalten
   @param ending     Formatierung: Endzeichen der Zeile  
*/
void smodel_Ak_print(FILE *file, smodel *smo, int k, char *tab,
		     char *separator, char *ending);

/**
   Schreiben der Gewichtsmatrix C eines Modells..
   @param file       Ausgabedatei
   @param smo        zu schreibendes Modell
   @param tab        zur Formatierung: fuehrende Tabs
   @param separator  Formatierung: Trennzeichen der Spalten
   @param ending     Formatierung: Endzeichen der Zeile  
*/
void smodel_C_print(FILE *file, smodel *smo, char *tab, char *separator, 
		    char *ending);
/**
   Schreiben der Mittelwertmatrix Mue eines Modells..
   @param file      Ausgabedatei
   @param smo       zu schreibendes Modell
   @param tab       zur Formatierung: fuehrende Tabs
   @param separator Formatierung: Trennzeichen der Spalten
   @param ending    Formatierung: Endzeichen der Zeile  
*/
void smodel_Mue_print(FILE *file, smodel *smo, char *tab, char *separator, 
		      char *ending);
/**
   Schreiben der Varianzenmatrix U eines Modells..
   @param file      Ausgabedatei
   @param smo       zu schreibendes Modell
   @param tab       zur Formatierung: fuehrende Tabs
   @param separator Formatierung: Trennzeichen der Spalten
   @param ending    Formatierung: Endzeichen der Zeile  
*/
void smodel_U_print(FILE *file, smodel *smo, char *tab, char *separator, 
			char *ending);
/**
   Schreiben des Anfangsbelegungsvektors eines Modells..
   @param file      Ausgabedatei
   @param smo       zu schreibendes Modell
   @param tab       zur Formatierung: fuehrende Tabs
   @param separator Formatierung: Trennzeichen der Spalten
   @param ending    Formatierung: Endzeichen der Zeile  
*/
void smodel_Pi_print(FILE *file, smodel *smo, char *tab, char *separator, 
		     char *ending);

void smodel_fix_print(FILE *file, smodel *smo, char *tab, char *separator, 
		     char *ending);
/**
   Schreiben der transponierten Uebergangsmatrix eines Modells..
   @param file      Ausgabedatei
   @param smo       zu schreibendes Modell
   @param tab       zur Formatierung: fuehrende Tabs
   @param separator Formatierung: Trennzeichen der Spalten
   @param ending    Formatierung: Endzeichen der Zeile  
*/
void smodel_Ak_print_transp(FILE *file, smodel *smo, int k, char *tab,
			    char *separator, char *ending);
/**
   Schreiben der transponierten Gewichtsmatrix eines Modells..
   @param file      Ausgabedatei
   @param smo       zu schreibendes Modell
   @param tab       zur Formatierung: fuehrende Tabs
   @param separator Formatierung: Trennzeichen der Spalten
   @param ending    Formatierung: Endzeichen der Zeile  
*/
void smodel_C_print_transp(FILE *file, smodel *smo, char *tab, char *separator, 
			   char *ending);
/**
   Schreiben der transponierten Mittelwertmatrix eines Modells..
   @param file   Ausgabedatei
   @param smo    zu schreibendes Modell
   @param tab    zur Formatierung: fuehrende Tabs
   @param separator  Formatierung: Trennzeichen der Spalten
   @param ending     Formatierung: Endzeichen der Zeile  
*/
void smodel_Mue_print_transp(FILE *file, smodel *smo, char *tab, 
			     char *separator, char *ending);
/**
   Schreiben der transponierten Varianzenmatrix eines Modells..
   @param file   Ausgabedatei
   @param smo    zu schreibendes Modell
   @param tab    zur Formatierung: fuehrende Tabs
   @param separator  Formatierung: Trennzeichen der Spalten
   @param ending     Formatierung: Endzeichen der Zeile  
*/
void smodel_U_print_transp(FILE *file, smodel *smo, char *tab, 
			       char *separator, char *ending);
/**
   Schreiben des transponierten Anfangsbelegungsvektors eines Modells..
   @param file   Ausgabedatei
   @param smo    zu schreibendes Modell
   @param tab    zur Formatierung: fuehrende Tabs
   @param separator  Formatierung: Trennzeichen der Spalten
   @param ending     Formatierung: Endzeichen der Zeile  
   */
void smodel_Pi_print_transp(FILE *file, smodel *smo, char *tab, char *ending);

void smodel_fix_print_transp(FILE *file, smodel *smo, char *tab, char *ending);

/**
   Schreiben eines SHMM-Modells in Matrixformat. 
   Input  smodel\_direct (heisst: Parameter als Matrix gespeichert.)
   @param file     Ausgabedatei
   @param smo_d    zu schreibendes Modell in Matrixformat
   @param multip   Modell wird in multip Kopien geschrieben
*/
void smodel_direct_print(FILE *file, smodel_direct *smo_d, int multip);

/** Alle Pointer des Modell freigeben und auf NULL setzen, Variablen 
    auf Null setzen..
    @param smo_d   SHMM-Modell-Struktur (\Ref{struct smodel_direct})
    @param check   Check-Struktur       (\Ref{struct shmm_check_t})
*/
void smodel_direct_clean(smodel_direct *smo_d, shmm_check_t *check); 

/** Teste Kompatibiliatet der Modell-Komponenten..
    @return ??
    @param smo_d  SHMM-Modell-Struktur  (\Ref{struct model_direct})
    @param check Check-Struktur         (\Ref{struct shmm_check_t})
*/
int smodel_direct_check_data(smodel_direct *smo_d, shmm_check_t *check); 

/** Berechnung c_im * b_im[omega]..
    @return       
    @param smo    continous HMM
    @param state  vorgegebener Zustand 
    @param m      Pointer auf Ergebnisvektor (vorher allozieren!)   
    @param omega  vorgegebener Ausgabewert 
*/
double smodel_calc_cmbm(smodel *smo, int state, int m, double omega);

/** Berechnung b[state][omega] (mixture-pdf)..
    @return       Wert
    @param smo    continous HMM
    @param state  vorgegebener Zustand
    @param omega  vorgegebener Ausgabewert
*/
double smodel_calc_b(smodel *smo, int state, double omega);

/** Computes probabilistic distance of two models
    @return the distance
    @param cm0  smodel used for generating random output
    @param cm   smodel to compare with
    @param maxT  maximum output length (for HMMs with absorbing states multiple
                 sequences with a toal langth of at least maxT will be 
		 generated)
    @param symmetric  flag, whether to symmetrize distance (not implemented yet)
    @param verbose  flag, whether to monitor distance in 40 steps. 
                    Prints to stdout (yuk!)
*/
double smodel_prob_distance(smodel *cm0, smodel *cm, int maxT, int symmetric, 
			    int verbose);

/** Berechnung c_im * B_im[omega] (Verteilungsfunktion!)..
    @return       
    @param smo    continous HMM
    @param state  vorgegebener Zustand 
    @param m      vorgegebener Index der Verteilung   
    @param omega  vorgegebener Ausgabewert
*/
double smodel_calc_cmBm(smodel *smo, int state, int m, double omega);

/** Berechnung B[state][omega] (mixture-Verteilungsfunktion!)..
    @return       
    @param smo    continous HMM
    @param state  vorgegebener Zustand
    @param omega  vorgegebener Ausgabewert
*/
double smodel_calc_B(smodel *smo, int state, double omega);

/** Intervall(a,b) mit  B(a) < 0.01, B(b) > 0.99
    @param smo    continous HMM
    @param state  vorgegebener Zustand
    @param a      return-value: left side
    @param b      return-value: right side
*/
void smodel_get_interval_B(smodel *smo, int state, double *a, double *b);


/** Verteilung der Daten und pro Zustand als Vektor von Daten und Gewichtungen..
    @return         0/-1
    @param smo      continous HMM
    @param sqd      sequence-field
    @param data     return: adress of the sorted data field
    @param weights  return: adress of the corresponding weight field
    @param len      return: vector-adress of the data lengths  
    @param wi       return: vector-adress of the state weights  
*/
int smodel_statedatacdf(smodel *smo, sequence_d_t *sqd, double ***data,
			double ***weights, int **len, double **wi);

/** Integraldifferenzen zwischen den Verteilungen der Daten und den theor.
    Verteilungsfunktionen pro Zustand.
    @return         0/-1
    @param smo      continous HMM
    @param sqd      sequence-field
    @param data     sorted data field          )
    @param weights  corresponding weight field ) = NULL -> internal calculation
    @param len      data lengths               )
    @param idiff    return: vector-adress of the distribution integral diffs 
*/
int smodel_idiff_index(smodel *smo, sequence_d_t *sqd, double **x, double **w,
		       int *len, double *wi, double **idiff);

/** Kolmogorov-Smirnov-Test ueber die Verteilungen der Daten und die 
    theor.Verteilungsfunktionen pro state.
    @return         0/-1
    @param smo      continous HMM
    @param sqd      sequence-field
    @param data     sorted data field            )
    @param weights  corresponding weight field   ) = NULL -> internal calcul.
    @param len      data lengths                 )
    @param wi       state weights                )
    @param prob     return: vector-adress of the KS-Test significance prob.
*/
int smodel_kstest_index(smodel *smo, sequence_d_t *sqd, double **x, double **w,
			int *len, double *wi, double **probks);

/** Chi-Square Test ueber die Verteilungen der Daten und die 
    theor.Verteilungsfunktionen pro state.
    @return         0/-1
    @param smo      continous HMM
    @param sqd      sequence-field
    @param data     sorted data field            )
    @param weights  corresponding weight field   ) = NULL -> internal calcul.
    @param len      data lengths                 )
    @param wi       state weights                )
    @param prob     return: vector-adress of the chi-square significance prob.s 
*/
int smodel_chisquare_index(smodel *smo, sequence_d_t *sqd, double **x, 
			   double **w, int *len, double *wi, double **prob);

int smodel_idiffsum_lowerbound(smodel *smo,sequence_d_t *sqd,double *idiffsum);

int smodel_count_free_parameter(smodel **smo, int smo_number);

#endif

/*@} (Doc++-Group: smodel) */
