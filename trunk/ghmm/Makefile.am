#  author       : Achim Gaedke
#  filename     : ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

SUBDIRS = ghmm tools tests doc

cvsclean: distclean
	$(RM) -f aclocal.m4 configure \
	Makefile.in stamp-h.in config.h.in \
	mkinstalldirs missing install-sh \
	doc/Makefile.in \
	ghmm/Makefile.in \
	tests/Makefile.in \
	tools/Makefile.in \
	$(distdir).tar.gz 

.PHONY: cvsclean

EXTRA_DIST = COPYRIGHT.txt GPL.txt LGPL.txt

COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

dist dist-all: copyright

copyright: distdir
	test -z "$(AWK)" && exit 0 ;\
	all_files=`cd $(top_distdir) && find . -type f` ; \
	for f in $$all_files ; do \
	  srcfile=$(top_srcdir)/$$f ; \
	  if grep '__copyright__ *$$' $$srcfile >/dev/null ; then \
	    distfile=$(top_distdir)/$$f ; \
	    rm $$distfile ; \
	    $(AWK) 'BEGIN {nr=0;while(getline < "$(COPYRIGHT)" ==1) copyright[nr++]=$$0;} /__copyright__ *$$/{sub("__copyright__.*","");for(i=0;i<nr;i++) print $$0 copyright[i];} //' \
	    $$srcfile >$$distfile; \
	  fi; \
	done
