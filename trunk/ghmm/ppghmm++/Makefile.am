#  author       : Peter Pipenbacher
#  filename     : $Source$
#  created      : 21.01.2002
#  $Id$
#
# __copyright__

# tests and examples

SUBDIRS = .
CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

INCLUDES = -I$(top_srcdir) $(XMLIO_INCLUDES)
LDADD    = $(top_builddir)/ppghmm++/libghmm++.a $(top_builddir)/ghmm/libghmm.a $(XMLIO_LIBS)

lib_LIBRARIES = libghmm++.a
libghmm___a_SOURCES = \
GHMM_DiscreteModel.cpp \
GHMM_ContinuousModel.cpp \
GHMM_AbstractModel.cpp \
GHMM_State.cpp \
GHMM_Sequences.cpp \
GHMM_ModelFactory.cpp \
GHMM_Document.cpp \
GHMM_Transition.cpp \
GHMM_State.cpp

pkginclude_HEADERS = \
GHMM_DiscreteModel.h \
GHMM_ContinuousModel.h \
GHMM_AbstractModel.h \
GHMM_State.h \
GHMM_Sequences.h \
GHMM_ModelFactory.h \
GHMM_Document.h \
GHMM_Transition.h \
GHMM_State.h


noinst_PROGRAMS = coin_toss_test chmm train_model generiere_nvs

coin_toss_test_SOURCES = coin_toss_test.cpp

chmm_SOURCES = chmm.cpp

train_model_SOURCES = train_model.cpp

generiere_nvs_SOURCES = generiere_nvs.cpp


COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

dist-hook:
	test -z "$(AWK)" && exit 0 ;\
	TMPFILE="$(distdir)/tmp"; \
	for f in $(DISTFILES) ; do \
	  distfile=$(distdir)/$$f ; \
	  if grep '__copyright__ *$$' $$distfile >/dev/null ; then \
	    mv $$distfile $$TMPFILE; \
	    $(AWK) 'BEGIN {nr=0;while(getline < "$(COPYRIGHT)" ==1) copyright[nr++]=$$0;} /__copyright__ *$$/{sub("__copyright__.*","");for(i=0;i<nr;i++) print $$0 copyright[i];} //' \
	    $$TMPFILE >$$distfile; \
	    touch -m -a -r $$TMPFILE $$distfile; \
	  fi; \
	done ; \
	rm -f $$TMPFILE
