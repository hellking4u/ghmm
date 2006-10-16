#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

SUBDIRS = .

SWIG_INTERFACE_FILES = ghmmwrapper.i wrapper_cseq.i wrapper_dseq.i wrapper_dpseq.i wrapper_cmodel.i wrapper_dmodel.i wrapper_xmlfile.i

all: setup.py ghmm.py ${SWIG_INTERFACE_FILES}
	python setup.py build

clean: setup.py
	python setup.py clean
	rm -f message.txt

install: setup.py all
	python setup.py install --prefix=${prefix}

distclean: setup.py
	python setup.py clean --all
	rm -f message.txt

.PHONY: all clean install distclean
