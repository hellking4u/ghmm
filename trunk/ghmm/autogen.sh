#! /bin/sh
#  author       : Achim Gaedke
#  filename     : ghmm/autogen.sh
#  created      : DATE: April 2001
#  $Id$

#move GNU m4 to head of PATH
for f in `type -ap m4` ; do
  echo `eval $f --version`|grep GNU >/dev/null && PATH=`dirname $f`:$PATH ;
done

#makes aclocal.m4 from acinclude.m4 and other files
aclocal

#scans configure.in and creates config.h.in
autoheader

#creates Makefile.in from Makefile.am
automake --add-missing

#creates configure from configure.in
autoconf


