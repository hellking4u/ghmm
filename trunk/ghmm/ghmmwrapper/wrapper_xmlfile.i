%{
#include <ghmm/xmlreader.h>
#include <ghmm/xmlwriter.h>
%}

typedef struct {

    int noModels;
    
    int modelType;
    
    union {
      ghmm_cmodel**  c;
      ghmm_dmodel**  d;
      ghmm_dpmodel** dp;
      ghmm_dsmodel** ds;
    } model;
    /* stupid swig */
} ghmm_xmlfile;

%extend ghmm_xmlfile {
        ghmm_cmodel*  get_cmodel(size_t index)  { return self->model.c[index]; }
        ghmm_dmodel*  get_dmodel(size_t index)  { return self->model.d[index]; }
        ghmm_dpmodel* get_dpmodel(size_t index) { return self->model.dp[index]; }
        ghmm_dsmodel* get_dsmodel(size_t index) { return self->model.ds[index]; }
}

extern ghmm_xmlfile* ghmm_xmlfile_parse(const char *filename);

extern int           ghmm_xmlfile_validate(const char *filename);

extern void          ghmm_xmlfile_write(ghmm_xmlfile* f, const char *file);

