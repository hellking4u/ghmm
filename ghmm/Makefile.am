#  author       : Achim Gaedke
#  filename     : ghmm/ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

MAINTAINERCLEANFILES = Makefile.in


INCLUDES = -I$(top_srcdir)
lib_LIBRARIES = libghmm.a

libghmm_a_SOURCES= foba.c foba.h \
                  matrix.c  matrix.h \
                  model.c model.h \
                  cluster.c cluster.h \
                  scluster.c scluster.h \
		  root_finder.c root_finder.h \
                  reestimate.c reestimate.h \
                  randvar.c randvar.h \
                  rng.c rng.h \
                  sequence.c sequence.h \
                  vector.c vector.h \
                  viterbi.c viterbi.h \
                  sfoba.c sfoba.h \
                  smodel.c smodel.h \
                  sreestimate.c sreestimate.h \
                  smixturehmm.c smixturehmm.h \
                  const.h \
                  mes.c mes.h \
                  mprintf.c mprintf.h \
                  scanner.c scanner.h \
                  gauss_tail.c gauss_tail.h \
                  smap_classify.c smap_classify.h

# excluded sgenerate, generation of sequences now in smodel
#                  sgenerate.c sgenerate.h \



pkginclude_HEADERS = foba.h \
                  matrix.h \
                  model.h \
                  cluster.h \
                  scluster.h \
		  root_finder.h \
                  reestimate.h \
                  randvar.h \
                  rng.h \
                  sequence.h \
                  vector.h \
                  viterbi.h \
                  sfoba.h \
                  smodel.h \
                  sreestimate.h \
                  smixturehmm.h \
                  const.h \
                  mes.h \
                  mprintf.h \
                  scanner.h \
                  gauss_tail.h \
                  smap_classify.h

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

