#  author       : Achim Gaedke
#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm++/examples/Makefile.am
#  created      : DATE: 28 June 2001
#  $Id$
#
# __copyright__

SUBDIRS= .
DISTCLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

INCLUDES = -I$(top_srcdir) $(XMLIO_INCLUDES)
LDADD = $(top_builddir)/ghmm++/libghmm++.a $(top_builddir)/ghmm/libghmm.a -lgxl $(XMLIO_LIBS)
EXTRADIST = ghmm.xml

noinst_PROGRAMS = sequence_test DiscretePD_test ghmm_test ghmm_convert
sequence_test_SOURCES = sequence_test.cpp
DiscretePD_test_SOURCES = DiscretePD_test.cpp
ghmm_test_SOURCES = ghmm_test.cpp
ghmm_convert_SOURCES = ghmm_convert.cpp

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





