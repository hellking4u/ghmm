#  author       : Achim Gaedke
#  filename     : ghmm/tests/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

include $(top_srcdir)/dist_hook.mak

CCLD     = $(CC)   -lstdc++

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
# in progress: test_sdfoba nullmodel ciscreen test_sdmodel 
noinst_PROGRAMS = randvar_test \
                  root_finder_test \
                  coin_toss_test \
                  two_states_three_symbols \
                  libxml-test \
                  chmm \
                  chmm_test \
                  shmm_viterbi_test \
                  test_gsl_ran_gaussian_tail \
                  sequences_test \
                  sequences_old_format \
                  label_higher_order_test \
		  read_fa

LDADD = $(top_builddir)/ghmm/.libs/libghmm.a
