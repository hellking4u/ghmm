#  author       : Achim Gaedke
#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm++/Makefile.am
#  created      : DATE: 28 June 2001
#  $Id$
#
# __copyright__

SUBDIRS= . examples
DISTCLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

INCLUDES = -I$(top_srcdir) $(XMLIO_INCLUDES)

pkgincludedir = $(includedir)/ghmm++

lib_LIBRARIES = libghmm++.a
libghmm___a_SOURCES = \
sequence.cpp \
sequences.cpp \
sequences_document.cpp \
InitialStates.cpp \
Emissions.cpp \
ghmm.cpp \
ghmm_graph.cpp \
ghmm_edge.cpp \
ghmm_alphabet.cpp \
ghmm_symbol.cpp

pkginclude_HEADERS= \
DiscretePD.h \
sequence.h \
sequences.h \
sequences_document.h \
InitialStates.h \
Emissions.h \
ghmm.h \
ghmm_graph.h \
ghmm_edge.h \
ghmm_alphabet.h \
ghmm_symbol.h

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





