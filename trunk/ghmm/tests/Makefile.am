#  author       : Achim Gaedke
#  filename     : /zpr/bspk/src/hmm/ghmm/tests/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

SUBDIRS = data

# tests and examples

AUTOMAKE_OPTIONS = dejagnu
EXPECT=$(EXPECT_FOR_DEJAGNU)
RUNTEST=$(RUNTEST_FOR_DEJAGNU)

CLEANFILES = ghmm.log ghmm.sum site.exp site.bak
MAINTAINERCLEANFILES = Makefile.in
# stuff from dejagnu

INCLUDES=-I$(top_srcdir)

#these tests will not be installed
noinst_PROGRAMS = root_finder_test coin_toss_test two_states_three_symbols chmm xml_parser_test shmm_viterbi_test generate_PHI test_gsl_ran_gaussian_tail sequences_test

LDADD = ../ghmm/libghmm.a

COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

# creates copyright tags and copies some files
dist-hook:
	DEJAGNU_FILES=`(cd $(srcdir) && find ghmm.tests lib config -name \*.exp)` ;\
	for f in $$DEJAGNU_FILES ; do \
	  dirpart=`dirname $$f`; \
	  test -d $(distdir)/$$dirpart || mkdir $(distdir)/$$dirpart ; \
	  cp -p $(srcdir)/$$f $(distdir)/$$dirpart ; \
	done; \
	test -z "$(AWK)" && exit 0 ;\
	TMPFILE="$(distdir)/tmp" ; \
	for f in $(DISTFILES) $$DEJAGNU_FILES ; do \
	  distfile=$(distdir)/$$f ; \
	  if grep '__copyright__ *$$' $$distfile >/dev/null ; then \
	    mv $$distfile $$TMPFILE; \
	    $(AWK) 'BEGIN {nr=0;while(getline < "$(COPYRIGHT)" ==1) copyright[nr++]=$$0;} /__copyright__ *$$/{sub("__copyright__.*","");for(i=0;i<nr;i++) print $$0 copyright[i];} //' \
	    $$TMPFILE >$$distfile; \
	    touch -m -a -r $$TMPFILE $$distfile; \
	  fi; \
	done ; \
	rm -f $$TMPFILE










