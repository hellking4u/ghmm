#ifndef SEQUENCE_H
#define SEQUENCE_H

/**@name Sequenzen (double und int) */
/*@{ (Doc++-Group: sequence) */
/** @name struct sequence_t
    Struktur zur Speicherung von Sequenzen. Eine Sequenz besteht aus
    einem int-Vektor und einem Sequenz Label. Sequenzen k\"onnen 
    unterschiedliche L\"ange haben.
 */
struct sequence_t {
  /** 
      Sequenzmatrix bestehend aus seq_number Zeilen. Die i-te Zeile
      hat die L\"ange seq_len[i].      
   */
  int **seq;
  /** Vektor der Sequenzl\"angen */
  int *seq_len;
  /** Vektor der Sequenzlabel */
  long *seq_label;
  /** Vektor der Sequenzid, double, da long evtl. nicht ausreichend */
  double *seq_id;
  /** Gewichtung der Squenzen; default 1 */
  double *seq_w;
  /** Anzahl der Sequenzen */
  long seq_number;
  /** Summe aller Sequenzgewichte, zur Sicherheit double */
  double total_w;
};
typedef struct sequence_t sequence_t;

/** @name struct sequence_d_t
    Struktur zur Speicherung von Sequenzen. Eine Sequenz besteht aus
    einem double Vektor und einem Sequenz Label. Sequenzen k\"onnen 
    unterschiedliche L\"ange haben
 */
struct sequence_d_t {
  /** 
      Sequenzmatrix bestehend aus seq_number Zeilen. Die i-te Zeile
      hat die L\"ange seq_len[i].      
   */
  double **seq;
  /** Vektor der Sequenzl\"angen */
  int *seq_len;
  /** Vektor der Sequenzlabel */
  long *seq_label;
  /** Vektor der Sequenzid, double, da long evtl. nicht ausreichend */
  double *seq_id;
  /** Gewichtung der Squenzen; default 1 */
  double *seq_w;
  /** Anzahl der Sequenzen */
  long seq_number;
  /** Summe aller Sequenzgewichte, zur Sicherheit double */
  double total_w;
};
typedef struct sequence_d_t sequence_d_t;

/*
  Wichtig: Include von (c)model.h zur Vermeidung von Fehlern beim Uebersetzen
  erst an dieser Stelle. (sequence.h und model.h includen sich gegenseitig)
*/

#include "model.h"
#include "smodel.h"

/**
   Einlesen von Sequenzen durch Aufruf von sequence\_read\_alloc..  
   @return pointer auf sequence_t
   @param filename    Filename
*/
sequence_t *sequence_read(char *filename);

/**
   Einlesen von und Speicheralloc fuer Sequenzen. Gibt verschiedene 
   Moeglichkeiten Sequenzen zu spezifizieren: direkte, explizite Angabe 
   (Keyword O) oder alle Sequenzen in lexikographischer Reihenfolge (Keyword L).
   @param p_seq_len  gefundene und dort allozierte Sequenzlaengen
   @param seq_number gefundene Sequenzanzahl
   @return ??
*/  
sequence_t *sequence_read_alloc(scanner_t *s);

/*
  Einlesen eines Sequenzfelds durch Aufruf von sequence\_d\_read\_alloc..  
  @return:           pointer auf sequence\_d\_t
  @param filename:   Filename
sequence_d_t *sequence_d_read_one(char *filename);
*/

/**
   Einlesen mehrerer Sequenzfelder durch wiederholten
   Aufruf von sequence\_d\_read\_alloc..  
   @return:           pointer auf Vektor von sequence\_d\_t
   @param filename:   Filename
   @param sqd_number: Anzahl der gelesenen Sequenzfelder
*/
sequence_d_t **sequence_d_read(char *filename, int *sqd_number);

/**
   Einlesen von und Speicheralloc fuer Double-Sequenzen. 
   Zur Zeit nur Explizite Angabe (Keyword O) realisiert.
   @return       eingelesene Sequenz (pointer auf sequence\_d\_t)
   @param s      Scanner-Struktur
*/  
sequence_d_t *sequence_d_read_alloc(scanner_t *s);

/** Truncate Sequences in a given sqd_field. useful for Testing;
   returns truncated sqd_field; 
   trunc_ratio 0: no truncation
   trunc_ratio 1: max. truncation
*/

sequence_d_t **sequence_d_truncate(sequence_d_t **sqd_in, int sqd_fields, 
				   double trunc_ratio, int seed);


/**
   Liefert lexikographisch geordnete Liste von Woertern der Laenge n
   aus einem Alphabet mit M Symbolen. In dieser Funktion ebenfalls alloc
   von sequence_t.
   @param n      Wortlaenge
   @param M      Alphabetgroesse
   @return Sequenz Struktur
*/
sequence_t *sequence_lexWords(int n, int M);

/**
   Bestimmt welches von vorgegebenen Modellen am besten zu einer vorgebenen
   Sequenz passt (hoechste Wahrscheinlichkeit). Liefert Index des 
   Modells zurueck und belegt log\_p mit der entsprechenden Wahrscheinlichkeit..
   @param mo            Vektor von Modellen
   @param model_number  Anzahl der Modelle in mo
   @param sequence      Sequenz
   @param seq_len       Sequenzlaenge
   @param log_p         log(Wahrscheinlichkeit, dass sequence von 
                        mo[model\_number] generiert wird)
   @return ??
*/
int sequence_best_model(model **mo, int model_number, int *sequence, 
			int seq_len, double *log_p);

/**
   Check, ob alle Symbole der Sequenzen erlaubt sind; d.h. aus [0 .. max\_symb - 1]..
   @param sq          zu pruefende Sequenzen
   @param max_symb    Anzahl der erlaubten Symbole; 
   @return            -1 (Fehler) / 0
*/
int sequence_check(sequence_t *sq, int max_symb);

/**
  Kopiert Symbole einer Sequenz; Speicher muss ausserhalb allokiert werden..
  @param target  Zielsequenz
  @param source  Quellsequenz
  @param len     Sequenzlaenge
  */
void sequence_copy(int *target, int *source, int len);

/**
  Kopiert Symbole einer Sequenz; Speicher muss ausserhalb allokiert werden..
  @param target  Zielsequenz
  @param source  Quellsequenz
  @param len     Sequenzlaenge
  */
void sequence_d_copy(double *target, double *source, int len);

/**
  Addiert alle Sequenzen aus source zu target hinzu. Allokiert Speicher 
  @param target  Zielsequenzen
  @param source  Quellsequenzen
  */
int sequence_add(sequence_t *target, sequence_t *source);

/**
  Addiert alle Sequenzen aus source zu target hinzu. Allokiert Speicher 
  @param target  Zielsequenzen
  @param source  Quellsequenzen
  */
int sequence_d_add(sequence_d_t *target, sequence_d_t *source);

/**
  Schreiben von Sequenzen..
  @param file        Ausgabedatei
  @param sequence    zu schreibende Sequenzen
  */
void sequence_print(FILE *file, sequence_t *sequence);

/** Schreiben von int Sequenzen in Mathematica Format (Liste von Listen).
 */
void sequence_mathematica_print(FILE *file, sequence_t *sq, char *name);

/**
  Schreiben von double-Sequenzen..
  @param file        Ausgabedatei
  @param seq_struct  Sequenzfeld
  @param discrete    0/1-Variable: 1 fuer int-Ausgabe
  */
void sequence_d_print(FILE *file, sequence_d_t *sqd, int discrete);

/** Schreiben von double Sequenzen in Mathematica Format (Liste von Listen).
 */
void sequence_d_mathematica_print(FILE *file, sequence_d_t *sqd, char *name);

/** GNUPLOT Format: ein Symbol pro Zeile, Seqs. durch Doppelleerzeile getrennt */
void sequence_d_gnu_print(FILE *file, sequence_d_t *sqd);
/**
  Nur free von Pointern innerhalb einer sequence_t struct,
  seq_number wird auf 0 gesetzt. Die Daten selbst bleiben erhalten!  
  @param seq_struct  structure, in der die Pointer freigegeben werden
  */
void sequence_clean(sequence_t *sq);

/**
  Nur free von Pointern innerhalb einer sequence_d_t struct,
  seq_number wird auf 0 gesetzt. Die Daten selbst bleiben erhalten!
  @param seq_struct  structure, in der die Pointer freigegeben werden
  */
void sequence_d_clean(sequence_d_t *sq);

/**
  free aller Pointer UND Speicherplaetze einer sequence_t struct.
  @param seq_struct  structure, in der Speicher freigegeben wird
  */
int sequence_free(sequence_t **sq);

/**
  free aller Pointer UND Speicherplaetze einer sequence_d_t struct.
  @param seq_struct  structure, in der Speicher freigegeben wird
  */
int sequence_d_free(sequence_d_t **sq);


/**
   liefert groesstes Symbol in einer Sequenz Struktur..
 */
int sequence_max_symbol(sequence_t *sq);

/**
   Alloc einer Sequenz Struktur fuer seq\_number Sequenzen. Es wird nur
   Speicher fuer die Pointer auf Sequenzen bereitgestellt, da die 
   Sequenzlaenge ja noch nicht bekannt ist. seq_number in sequence_t wird
   richtig belegt. Die Vektoren seq\_len und seq\_label werden auch allociert
   und mit 0 initialisiert.
   @param seq_number:  Anzahl der Sequenzen
   @return:            leere Sequenz Struktur.
*/
sequence_t *sequence_calloc(long seq_number);

/**
   Gleiche Funktion wie sequence_calloc fuer sequence_d_t
 */
sequence_d_t *sequence_d_calloc(long seq_number);

/**
   Erzeugen einer Double-Kopie eines Int-Sequenzfeldes..
   @return        zu erzeugendes double-Sequenzfeld
   @param sq      int-Sequenzfeld, das kopiert wird
   */
sequence_d_t *sequence_d_create_from_sq(const sequence_t *sq);

/**
   Erzeugen einer Int-Kopie eines Double-Sequenzfeldes..
   @return        zu erzeugendes int-Sequenzfeld
   @param sqd     double-Sequenzfeld, das kopiert wird
   */
sequence_t *sequence_create_from_sqd(const sequence_d_t *sqd);

/** 
  Maximale Sequenzl\"ange eines Sequenzfeldes ermitteln..
 */
int sequence_d_max_len(const sequence_d_t *sqd);

/**
  Berechnung der Mittelwertsequenz. Falls Sequenzen unterschiedlicher
  L\"ange vorliegen, wird die gr\"o\ss{}te L\"ange ermittelt und die
  fehlenden Teile der k\"urzeren Sequenzen werden mit null aufgef\"ullt.
  */
sequence_d_t *sequence_d_mean(const sequence_d_t *sqd);

/**
  Berechnung der Scattermatrix. Falls Sequenzen unterschiedlicher
  L\"ange vorliegen, wird die gr\"o\ss{}te L\"ange ermittelt und die
  fehlenden Teile der k\"urzeren Sequenzen werden NICHT BERUECKSICHTIGT!.
  @return        zu erzeugende double-Matrix
  @param sqd     double-Sequenzfeld
  @param sqd     Gr\"o\ss{}e der quadratischen Matrix
  */
double **sequence_d_scatter_matrix(const sequence_d_t *sqd, int *dim);

/**
  Schreiben von double-Sequenzen als Input-File fuer das K-Means-Programm
  @param filename        Ausgabedatei
  @param sqd             Sequenzfeld
  @param cluster_number  Anzahl der zu bildenden Cluster
  @param start_partition 1: Labels als Start-Partition der Clusterung nehmen
*/
int sequence_d_kmeans_print(char *filename, sequence_d_t *sqd, 
			    int cluster_number, int start_partition);

/**
  Lesen des Ergebnisses einer kmeans-Clusterung, d.h. zugeordneten Cluster werden 
  als Labels der Sequenzen gesetzt..
  @param filename        Ausgabedatei
  @param sqd             Sequenzfeld
  @param max_label       return: groesstes gelesenes Label
*/
int sequence_d_read_kmeans_label(char *filename, sequence_d_t *sqd, 
				 int *max_label);

/** liefert den ersten Index der Sequenz, deren ID mit der uebergebenen ID
    uebereinstimmt */
long sequence_d_index(sequence_d_t *sqd, double id);

/** liefert den ersten Index der Sequenz, deren ID mit der uebergebenen ID
    uebereinstimmt */
long sequence_index(sequence_t *sq, double id);

/**
   Funktion zur Klasseneinteilung einer Sequenz. Spar- und Tilgungsphase
   sind getrennt, Sondersymbole werden nicht summiert.
   Details auch im Source Kommentar.
   @param O Sequenz
   @param index Sequenzindex
 */
int sequence_d_class(const double *O, int index, double *osum, int *tilgphase);

/** Ist O[index] eine SPE Ausgabe?. Unterscheidung Tilg. und SPE
   ueber vorangegangenes symbol_zut */
int sequence_d_is_spar(double *O, int index);


/** Ist O[index] eine Tilg. Ausgabe? */

int sequence_d_is_tilg(double *O, int index);


/** zufaellige Aufteilung eines Seq. Feldes in Test und Trainingssequenzen */
int sequence_d_partition(sequence_d_t *sqd, sequence_d_t * sqd_train, 
			 sequence_d_t *sqd_test, double train_ratio);


/** kopiert source[s_num] nach target[t_num] inkl. label, id, w, len.
    Keine allocs hier!
*/
void sequence_d_copy_all(sequence_d_t *target, long t_num, 
			 sequence_d_t *source, long s_num);

int sequence_d_labels_from_kmeans(sequence_d_t *sqd, int smo_number);


int sequence_d_is_sonder(double O);

/** Likelihood function in a mixture model:
   $sum_k w^k log( sum_c (alpha_c p(O^k | lambda_c)))$
*/

int sequence_d_mix_like(smodel **smo, int  smo_number, sequence_d_t *sqd, double *like);

int sequence_d_mix_likeBS(smodel **smo, int  smo_number, sequence_d_t *sqd,
			  double *like);

#endif
/*@} (Doc++-Group: sequence) */
