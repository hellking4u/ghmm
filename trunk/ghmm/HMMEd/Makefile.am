#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id: Makefile.am 1828 2007-05-16 12:27:49Z grunau $
#
# __copyright__

SUBDIRS = .

PYTHON=python

PYTHON_FILES = HMMEd \
               HMMEditor.py \
               ObjectHMM.py

all: setup.py ${PYTHON_FILES}
	$(PYTHON) setup.py build

clean: setup.py
	$(PYTHON) setup.py clean

install: setup.py all
	$(PYTHON) setup.py install --prefix=${prefix}

distclean: setup.py
	$(PYTHON) setup.py clean --all

.PHONY: all clean install distclean
