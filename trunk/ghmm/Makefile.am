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

EXTRA_DIST = COPYRIGHT.txt GPL.txt LGPL.txt ToDo

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
