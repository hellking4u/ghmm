#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

all: setup.py
	swig -noruntime -python -nodefault ghmmwrapper.i
	python setup.py build

clean: setup.py
	python setup.py clean

install: setup.py
	swig -noruntime -python -nodefault ghmmwrapper.i
	python setup.py build
	python setup.py install --prefix=${prefix}


.PHONY: all clean install 
