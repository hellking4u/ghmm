#ghmm/Makefile.am

SUBDIRS = ghmm tools tests doc

cvsclean: distclean
	$(RM) -f aclocal.m4 configure \
	Makefile.in stamp-h.in config.h.in \
	mkinstalldirs missing install-sh \
	doc/Makefile.in \
	ghmm/Makefile.in \
	tests/Makefile.in \
	tools/Makefile.in

.PHONY: cvsclean

EXTRA_DIST = COPYRIGHT.txt GPL.txt LGPL.txt

COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

dist dist-all: copyright

copyright: distdir
	cr_files="`(cd $(top_distdir) && find . -name \*.h -o -name \*.c)` tools/ghmm-config.in" ; \
	for f in $$cr_files ;do \
	  distfile=$(top_distdir)/$$f ; \
	  srcfile=$(top_srcdir)/$$f ; \
	  rm $$distfile ; \
	  awk '/__copyright__/{system("cat $(COPYRIGHT)");next} //' $$srcfile >$$distfile; \
	done
