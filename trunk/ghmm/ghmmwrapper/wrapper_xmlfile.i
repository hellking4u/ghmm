%{
#include <ghmm/xmlreader.h>
#include <ghmm/xmlwriter.h>
%}

struct ghmm_xmlfile {

    int noModels;
    
    int modelType;
    
    union {
      ghmm_cmodel * * c;
      ghmm_dmodel * * d;
      ghmm_dpmodel * * dp;
      ghmm_dsmodel * * ds;
    } model;
    /* stupid swig */
};
typedef struct ghmm_xmlfile ghmm_xmlfile;

extern ghmm_xmlfile* ghmm_xmlfile_parse(const char *filename);

extern int           ghmm_xmlfile_validate(const char *filename);

extern void          ghmm_xmlfile_write(ghmm_xmlfile* f, const char *file);

