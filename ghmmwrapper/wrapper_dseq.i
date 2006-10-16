/*==========================================================================
  ===== discrete sequences ================================================= */
typedef struct ghmm_dseq {
  /** sequence array. sequence[i] [j] = j-th symbol of i-th seq.
   */
  int **seq;

  /** matrix of state ids, can be used to save the viterbi path during sequence generation.
      ATTENTION: is NOT allocated by ghmm_dseq_calloc  */
  int **states;

  /** array of sequence length */
  int *seq_len;
#ifdef GHMM_OBSOLETE
  /**  array of sequence labels */
  long *seq_label;
#endif /* GHMM_OBSOLETE */
  /**  array of sequence IDs*/
  double *seq_id;
  /** positiv! sequence weights.  default is 1 = no weight */
  double *seq_w;
  /** total number of sequences */
  long seq_number;
  /** sum of sequence weights */
  double total_w;

  /* matrix of state labels corresponding to seq */
  int **state_labels;
  /* number of labels for each sequence */
  int *state_labels_len;
} ghmm_dseq;

extern int ghmm_dseq_free(ghmm_dseq **sq);
extern ghmm_dseq* ghmm_dseq_calloc(long number);

%delobject ghmm_dseq::subseq_free;
%delobject ghmm_dseq_subseq_free;

%extend ghmm_dseq {
        ghmm_dseq(long number) { return ghmm_dseq_calloc(number); }
        ~ghmm_dseq() { ghmm_dseq_free(&self); }

        int calloc_state_labels(); 

        ghmm_dseq* get_singlesequence(int index);

        int subseq_free();

        int max_symbol();

        int add(ghmm_dseq *source);

        int check(int max_symb);

        void clean();

        int* getSequence(int index) { return self->seq[index]; }
        void setSequence(int seqno, int *O) { self->seq[seqno] = O; }
        int getSymbol(int seqno, int index) { return self->seq[seqno][index]; }
        void setSymbol(int seqno, int index, int value) { self->seq[seqno][index] = value; }

        int  getLength(int i) { return self->seq_len[i]; }
        void setLength(int i, int len) { self->seq_len[i] = len; }

        int  getLabelsLength(int i) { return self->state_labels_len[i]; }
        void setLabelsLength(int i, int len) { self->state_labels_len[i] = len; }

        double getWeight(int i) { return self->seq_w[i]; }
        void   setWeight(int i, double w) { self->seq_w[i] = w; }

        void copyStateLabel(int index, ghmm_dseq *target, int no)
            {
                int length = self->state_labels_len[index];
                target->state_labels_len[no] = length;
                target->state_labels[no]= malloc(self->state_labels_len[index] * sizeof(int));
                memcpy(target->state_labels[no], self->state_labels[index], length);
            }

        void write(char* filename)
            {
                FILE* file;
                file = fopen(filename, "at");
                ghmm_dseq_print(self, file);
                fclose(file);
            }
}

REFERENCE_ARRAY(ghmm_dseq, dseq_ptr)