#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

SUBDIRS = .
MAINTAINERCLEANFILES = Makefile.in

#NOTE: C-Files will be seperatly compiled by Python setup
DISTFILES = ghmm.py ghmmhelper.py \
            ghmmwrapper.i \
            setup.py \
            sdclass_change.c sdclass_change.h read_cxml.c read_cxml.h gql.c \
            mixture.py \
            test100.sqd test2.smo


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

