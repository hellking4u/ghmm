#  author       : Peter Pipenbacher
#  filename     : $Source$
#  created      : 21.01.2002
#  $Id$
#
#  revision date   : $Date$
#  ___copyright__
# 

# tests and examples

SUBDIRS = . examples
CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

INCLUDES = -I$(top_srcdir) $(XMLIO_INCLUDES)

pkgincludedir = $(includedir)/ghmm++

LDFLAGS = $(XMLIO_LIBS)

lib_LTLIBRARIES = libghmm++.la
libghmm___la_LDFLAGS = -version-info @GHMM_VERSION_CURRENT@:@GHMM_VERSION_REVISION@:@GHMM_VERSION_AGE@
libghmm___la_SOURCES = \
GHMM_AbstractModel.cpp \
GHMM_DiscreteModel.cpp \
GHMM_ContinuousModel.cpp \
GHMM_State.cpp \
GHMM_Sequence.cpp \
GHMM_Sequences.cpp \
GHMM_Document.cpp \
GHMM_Transition.cpp \
GHMM_Emission.cpp \
GHMM_IntVector.cpp \
GHMM_DoubleVector.cpp \
GHMM_DoubleMatrix.cpp \
GHMM_IntMatrix.cpp \
GHMM_Alphabet.cpp \
GHMM_Toolkit.cpp \
GHMM_GMLDiscreteModel.cpp \
GHMM_GMLContinuousModel.cpp \
GHMM_GMLState.cpp \
GHMM_GMLTransition.cpp \
GHMM_GMLAlphabet.cpp \
GHMM_GMLDoc.cpp \
GHMM_convertXMLtoC.cpp 

pkginclude_HEADERS = \
begin_code.h \
close_code.h \
GHMM_DiscreteModel.h \
GHMM_ContinuousModel.h \
GHMM_AbstractModel.h \
GHMM_State.h \
GHMM_Sequence.h \
GHMM_Sequences.h \
GHMM_Document.h \
GHMM_Transition.h \
GHMM_Emission.h \
GHMM_State.h \
GHMM_Types.h \
GHMM_IntVector.h \
GHMM_DoubleVector.h \
GHMM_DoubleMatrix.h \
GHMM_IntMatrix.h \
GHMM_Toolkit.h \
GHMM_Alphabet.h \
GHMM.h \
GHMM_convertXMLtoC.h \
GHMM_GMLDiscreteModel.h \
GHMM_GMLContinuousModel.h \
GHMM_GMLState.h \
GHMM_GMLTransition.h \
GHMM_GMLAlphabet.h



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
