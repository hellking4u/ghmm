#*******************************************************************************
#
#       This file is part of the General Hidden Markov Model Library,
#       GHMM version __VERSION__, see http://ghmm.org
#
#       Filename: ghmm/ghmm/Makefile.am
#       Authors:  Wasinee Rungsarityotin, Benjamin Georgi, Alexander Schliep,
#                 Benjamin Rich, Janne Grunau
#
#       Copyright (C) 1998-2004 Alexander Schliep
#       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
#       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
#                               Berlin
#
#       Contact: schliep@ghmm.org
#
#       This library is free software; you can redistribute it and/or
#       modify it under the terms of the GNU Library General Public
#       License as published by the Free Software Foundation; either
#       version 2 of the License, or (at your option) any later version.
#
#       This library is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#       Library General Public License for more details.
#
#       You should have received a copy of the GNU Library General Public
#       License along with this library; if not, write to the Free
#       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#
#       This file is version $Revision$
#                       from $Date$
#             last change by $Author$.
#
#******************************************************************************/

SUBDIRS = .
MAINTAINERCLEANFILES = Makefile.in

INCLUDES = -I$(top_srcdir) @GHMM_CFLAGS@
lib_LTLIBRARIES = libghmm.la

libghmm_la_LDFLAGS = -version-info @GHMM_VERSION_CURRENT@:@GHMM_VERSION_REVISION@:@GHMM_VERSION_AGE@
libghmm_la_SOURCES= ghmm.h \
                    rng.c rng.h \
                    randvar.c randvar.h \
                    gauss_tail.c \
                    mes.c mes.h \
                    ghmm_internals.c ghmm_internals.h \
                    mprintf.c mprintf.h \
                    vector.c vector.h \
                    matrix.c  matrix.h \
                    scanner.c scanner.h \
                    linkedlist.c \
                    sequence.c sequence.h \
                    psequence.c psequence.h \
                    xmlreader.c xmlreader.h \
                    model.c model.h \
                    foba.c foba.h \
                    viterbi.c viterbi.h \
                    reestimate.c reestimate.h \
                    gradescent.c gradescent.h \
                    kbest.c kbest.h \
                    discrime.c discrime.h \
                    modelutil.c \
                    cluster.c cluster.h \
                    root_finder.c \
                    smodel.c smodel.h \
                    sfoba.c sfoba.h \
                    sviterbi.c sviterbi.h \
                    sreestimate.c sreestimate.h \
                    scluster.c scluster.h \
                    sgenerate.c sgenerate.h \
                    smixturehmm.c smixturehmm.h \
                    smap_classify.c smap_classify.h \
                    sdmodel.c sdmodel.h \
                    sdfoba.c sdfoba.h \
                    sdviterbi.c \
                    pmodel.c pmodel.h \
                    pviterbi.c pviterbi.h \
                    pviterbi_propagate.c pviterbi_propagate.h


pkginclude_HEADERS = ghmm.h \
                  sequence.h \
		  psequence.h \
                  xmlreader.h \
                  model.h \
		  foba.h \
                  viterbi.h \
                  gradescent.h \
                  kbest.h \
                  discrime.h \
                  cluster.h \
                  scluster.h \
                  reestimate.h \
                  sfoba.h \
                  sviterbi.h \
                  smodel.h \
		  sdmodel.h \
                  sdfoba.h \
                  sreestimate.h \
                  sgenerate.h \
                  smixturehmm.h \
                  smap_classify.h \
		  pmodel.h \
		  pviterbi.h \
		  pviterbi_propagate.h \
		  rng.h \
		  scanner.h \
		  obsolete.h \
		  unsupported.h

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

