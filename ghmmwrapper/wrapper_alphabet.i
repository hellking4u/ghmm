typedef struct {
  int id;
  char* description;
  unsigned int size;
  char** symbols;
} ghmm_alphabet;

%extend ghmm_alphabet {
        ghmm_alphabet(size_t size, char* desc) {
            ghmm_alphabet* a = calloc(1, sizeof(ghmm_alphabet));
            a->symbols = calloc(size, sizeof(char*));
            a->description = desc;
            return a;
        }
       ~ghmm_alphabet() {
            int i;
            for(i=0; i<self->size; ++i)
                free(self->symbols[i]);
            free(self->symbols);
            free(self->description);
            free(self);
        }

        char* getSymbol(size_t index) { return self->symbols[index]; }
        void  setSymbol(char* s, size_t index) { self->symbols[index] = s; }
}