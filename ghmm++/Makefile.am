#  author       : Peter Pipenbacher
#  filename     : $Source$
#  created      : 21.01.2002
#  $Id$
#
# __copyright__

# tests and examples

SUBDIRS = . examples
CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

INCLUDES = -I$(top_srcdir) $(XMLIO_INCLUDES)

pkgincludedir = $(includedir)/ghmm++

lib_LIBRARIES = libghmm++.a
#libghmm___la_LDFLAGS = -version-info @GHMM_VERSION_CURRENT@:@GHMM_VERSION_REVISION@:@GHMM_VERSION_AGE@
libghmm___a_SOURCES = \
GHMM_DiscreteModel.cpp \
GHMM_ContinuousModel.cpp \
GHMM_AbstractModel.cpp \
GHMM_State.cpp \
GHMM_Sequence.cpp \
GHMM_Sequences.cpp \
GHMM_ModelFactory.cpp \
GHMM_Document.cpp \
GHMM_Transition.cpp \
GHMM_Emission.cpp \
GHMM_State.cpp \
GHMM_IntVector.cpp \
GHMM_DoubleVector.cpp \
GHMM_DoubleMatrix.cpp

pkginclude_HEADERS = \
begin_code.h \
close_code.h \
GHMM_DiscreteModel.h \
GHMM_ContinuousModel.h \
GHMM_AbstractModel.h \
GHMM_State.h \
GHMM_Sequence.h \
GHMM_Sequences.h \
GHMM_ModelFactory.h \
GHMM_Document.h \
GHMM_Transition.h \
GHMM_Emission.h \
GHMM_State.h \
GHMM_Types.h \
GHMM_IntVector.h \
GHMM_DoubleVector.h \
GHMM_DoubleMatrix.h


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
