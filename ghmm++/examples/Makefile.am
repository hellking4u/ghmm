#  author       : Achim Gaedke
#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm++/examples/Makefile.am
#  created      : DATE: 28 June 2001
#  $Id$
#  revision date : $Date$
# 
# __copyright__

SUBDIRS= .
CLEANFILES = \
chmm.tmp \
message.txt
DISTCLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

#
# if you want to link ghmm++ from C-code, you must use dynamic linking, 
#
# NOTICE: It is important to link libstdc++ when linking from gcc to a C++ library. 
#
CCLD     = $(CC)   -lstdc++
CXXLD    = $(CXX)  -lstdc++

INCLUDES = -I$(top_srcdir) $(XMLIO_INCLUDES)
LDADD    = $(top_builddir)/ghmm++/libghmm++.la $(top_builddir)/ghmm/libghmm.la
EXTRA_DIST = nvs.xml nvs.mod nvs.sqd gene_hmm.gml

noinst_PROGRAMS = dna chmm generiere_nvs generiere_nvs_old_files coin_toss_test vector_test create_sequence read_graphml
 
coin_toss_test_SOURCES = coin_toss_test.cpp
chmm_SOURCES = chmm.cpp
generiere_nvs_SOURCES = generiere_nvs.cpp
generiere_nvs_old_files_SOURCES = generiere_nvs_old_files.cpp
dna_SOURCES = dna.cpp
vector_test_SOURCES = vector_test.cpp
create_sequence_SOURCES = create_sequence.cpp
read_graphml_SOURCES = read_graphml.cpp
# read_cxml_SOURCES = read_cxml.c


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





