#  author       : Achim Gaedke
#  filename     : ghmm/tests/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

# tests and examples

AUTOMAKE_OPTIONS = dejagnu

CLEANFILES = ghmm.log ghmm.sum site.exp site.bak
MAINTAINERCLEANFILES = Makefile.in

INCLUDES=-I$(top_srcdir)

#these tests will not be installed
noinst_PROGRAMS = root_finder_test coin_toss_test two_states_three_symbols

LDADD = ../ghmm/libghmm.a

COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

dist-hook:
	test -z "$(AWK)" && exit 0 ;\
	for f in $(DISTFILES) ; do \
	  distfile=$(distdir)/$$f ; \
	  if grep '__copyright__ *$$' $$distfile >/dev/null ; then \
	    mv $$distfile tmp; \
	    $(AWK) 'BEGIN {nr=0;while(getline < "$(COPYRIGHT)" ==1) copyright[nr++]=$$0;} /__copyright__ *$$/{sub("__copyright__.*","");for(i=0;i<nr;i++) print $$0 copyright[i];} //' \
	    tmp >$$distfile; \
	    touch -m -a -r tmp $$distfile; \
	  fi; \
	done ; \
	rm -f tmp




