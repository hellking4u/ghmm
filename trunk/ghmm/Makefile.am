#  author       : Achim Gaedke
#  filename     : ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

SUBDIRS = @SUBDIRS@

MAINTAINERCLEANFILES = \
aclocal.m4 \
configure \
Makefile.in \
stamp-h.in \
config.h.in \
mkinstalldirs \
missing \
install-sh \
INSTALL \
COPYING \
$(distdir).tar.gz -r $(distdir) \
ltmain.sh \
config.guess \
config.sub \
$(distdir).tar.gz \
-r $(distdir)

EXTRA_DIST = \
COPYRIGHT.txt \
GPL.txt \
LGPL.txt \
ToDo \
win_config.h \
README.Win

COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

dist-hook:
	test -z "$(AWK)" && exit 0 ;\
	TMPFILE="$(distdir)/tmp" ; \
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
