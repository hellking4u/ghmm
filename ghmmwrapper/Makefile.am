#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

SUBDIRS = .

PYTHON=python

PYTHONVER := $(shell $(PYTHON) -V 2>&1 | cut -c 8-10)

SWIG_INTERFACE_FILES = ghmmwrapper.i \
                       wrapper_alphabet.i \
                       wrapper_cseq.i \
                       wrapper_dseq.i \
                       wrapper_dpseq.i \
                       wrapper_cmodel.i \
                       wrapper_dmodel.i \
                       wrapper_dpmodel.i \
                       wrapper_xmlfile.i

all: setup.py ghmm.py ${SWIG_INTERFACE_FILES}
	$(PYTHON) setup.py build

clean: setup.py
	$(PYTHON) setup.py clean
	rm -f message.txt

install: setup.py all
	$(PYTHON) setup.py install --prefix=${prefix}

test: ghmmunittests.py all
	PYTHONPATH=./build/lib.linux-i686-$(PYTHONVER)/:$(PYTHONPATH) LD_LIBRARY_PATH=../ghmm/.libs:$(LD_LIBRARY_PATH) $(PYTHON) ghmmunittests.py

coverage: ghmmunittests.py all
	coverage.py -e
	coverage.py -x ghmmunittests.py
	coverage.py -r ghmm.py

distclean: setup.py
	$(PYTHON) setup.py clean --all
	rm -f message.txt

.PHONY: all clean install distclean test coverage
