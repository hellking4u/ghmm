#  author       : Achim Gaedke
#  filename     : ghmm/tools/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

MAINTAINERCLEANFILES = Makefile.in

BUILT_SOURCES = 
INCLUDES = -I$(top_srcdir)
bin_PROGRAMS = probdist cluster scluster

probdist_SORUCES = probdist.c
cluster_SOURCES  = cluster.c
scluster_SOURCES = scluster.c

LDADD = ../ghmm/libghmm.a
bin_SCRIPTS = ghmm-config

COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

dist-hook:
	test -z "$(AWK)" && exit 0 ;\
	TMPFILE=`mktemp $(distdir)/dist-hook.XXXXXX`; \
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






