#! /bin/sh
#  author       : Achim Gaedke
#  filename     : ghmm/autogen.sh
#  created      : DATE: April 2001
#  $Id$


#makes aclocal.m4 from acinclude.m4 and other files
aclocal

#scans configure.in and creates config.h.in
autoheader

#creates Makefile.in from Makefile.am
automake --add-missing --copy

#creates configure from configure.in
autoconf
