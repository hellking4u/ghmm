#  author       : Achim Gaedke
#  filename     : ghmm/ghmm++/Makefile.am
#  created      : DATE: 28 June 2001
#  $Id$
#
# __copyright__

SUBDIRS=
DISTCLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

INCLUDES = -I$(top_srcdir)
LDADD = libghmm++.a $(top_builddir)/ghmm/libghmm.a -lxmlio -lexpat

bin_PROGRAMS = sequences_reader sequence_test
sequences_reader_SOURCES = sequences_reader.cpp
sequence_test_SOURCES = sequence_test.cpp

lib_LIBRARIES = libghmm++.a
libghmm___a_SOURCES = \
sequences.cpp sequences.h \
XMLIO_ArrayObject.cpp XMLIO_ArrayObject.h

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





