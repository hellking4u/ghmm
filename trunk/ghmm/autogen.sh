#! /bin/sh

#clean up
rm -f aclocal.m4 configure \
Makefile.in stamp-h.in config.h.in \
mkinstalldirs missing install-sh \
doc/Makefile.in \
ghmm/Makefile.in \
tests/Makefile.in \
tools/Makefile.in

#makes aclocal.m4 from acinclude.m4 and other files
aclocal

#scans configure.in and creates config.h.in
autoheader

#creates Makefile.in from Makefile.am
automake --add-missing --copy

#creates configure from configure.in
autoconf
