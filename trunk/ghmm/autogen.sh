#!/bin/sh

#clean up
rm -f aclocal.m4 configure Makefile.in Makefile stamp-h.in stamp-h \
mkinstalldirs missing install-sh \
config.h config.h.in \
config.cache config.log config.status \
doc/Makefile doc/Makefile.in doc/docu.dxx \
ghmm/Makefile ghmm/Makefile.in \
tests/Makefile tests/Makefile.in \
tools/Makefile tools/Makefile

#makes aclocal.m4 from acinclude.m4 and other files
aclocal

#scans configure.in and creates config.h.in
autoheader

#creates Makefile.in from Makefile.am
automake --add-missing --copy

#creates configure from configure.in
autoconf
