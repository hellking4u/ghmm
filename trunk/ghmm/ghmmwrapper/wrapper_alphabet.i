typedef struct {
  int id;
  char* description;
  unsigned int size;
  char** symbols;
} ghmm_alphabet;

%extend ghmm_alphabet {
        char* getSymbol(size_t index) { return self->symbols[index]; }
}