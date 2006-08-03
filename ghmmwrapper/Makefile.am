#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

SUBDIRS = .

all: setup.py ghmmwrapper.i ghmm.py
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
