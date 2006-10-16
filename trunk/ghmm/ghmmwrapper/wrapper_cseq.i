/*==========================================================================
  ===== continous sequences ================================================ */

typedef struct ghmm_cseq {
  /** sequence array. sequence[i][j] = j-th symbol of i-th seq. */
  double **seq;
  /** array of sequence length */
  int *seq_len;
#ifdef GHMM_OBSOLETE
  /**  array of sequence labels */
    long *seq_label;
#endif /* GHMM_OBSOLETE */
  /**  array of sequence IDs*/
  double *seq_id;
  /** positive! sequence weights.  default is 1 = no weight */
  double *seq_w;
  /** total number of sequences */
  long seq_number;
  /** sum of sequence weights */
  double total_w;
} ghmm_cseq;

extern int ghmm_cseq_free(ghmm_cseq **csq);
extern ghmm_cseq* ghmm_cseq_calloc(long number);

%delobject ghmm_cseq::subseq_free;
%delobject ghmm_cseq_subseq_free;

%extend ghmm_cseq {
        ghmm_cseq(long number) { return ghmm_cseq_calloc(number); }
        ~ghmm_cseq() { ghmm_cseq_free(&self); }

        //ghmm_cseq** truncate(int sqd_fields, double trunc_ratio, int seed);

        ghmm_cseq* get_singlesequence(int index);

        int subseq_free();

        int add(ghmm_cseq *source);

        void clean();

        int partition(ghmm_cseq *sqd_train, ghmm_cseq *sqd_test, double train_ratio);

        void copy_all(long t_num, ghmm_cseq *source, long s_num);

        double* getSequence(int index) { return self->seq[index]; }
        void setSequence(int seqno, double *O) { self->seq[seqno] = O; }
        double getSymbol(int seqno, int index) { return self->seq[seqno][index]; }
        void setSymbol(int seqno, int index, double value) { self->seq[seqno][index] = value; }

        int  getLength(int i) { return self->seq_len[i]; }
        void setLength(int i, int len) { self->seq_len[i] = len; }

        double getWeight(int i) { return self->seq_w[i]; }
        void   setWeight(int i, double w) { self->seq_w[i] = w; }

        void write(char* filename, int discrete)
            {
                FILE* file;
                file = fopen(filename, "at");
                ghmm_cseq_print(self, file, discrete);
                fclose(file);
            }
}