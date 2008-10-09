#  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id$
#
# __copyright__

include $(top_srcdir)/dist_hook.mak

SUBDIRS = .

SWIG_INTERFACE_FILES = ghmmwrapper.i \
                       wrapper_alphabet.i \
		       wrapper_arrays.i \
                       wrapper_cseq.i \
                       wrapper_dseq.i \
                       wrapper_dpseq.i \
                       wrapper_cmodel.i \
                       wrapper_dmodel.i \
                       wrapper_dpmodel.i \
                       wrapper_xmlfile.i

EXTENSION_C_FILES = gql.c \
                    pclasschange.c \
                    pclasschange.h \
                    sclass_change.c \
                    sclass_change.h \
                    sdclass_change.c \
                    sdclass_change.h

PYTHON_FILES = class_change.py \
               ghmmhelper.py \
               ghmm.py \
               ghmmunittests.py \
               mixture.py \
               modhmmer.py \
               setup.py

PYTHON_GATO_FILES = ghmm_gato/DataStructures.py \
                    ghmm_gato/EditObjectAttributesDialog.py \
                    ghmm_gato/GatoGlobals.py \
                    ghmm_gato/Graph.py \
                    ghmm_gato/GraphUtil.py \
                    ghmm_gato/HMMEd.py \
                    ghmm_gato/__init__.py \
                    ghmm_gato/xmlutil.py

EXTRA_DIST = ${SWIG_INTERFACE_FILES} ${EXTENSION_C_FILES} ${PYTHON_FILES} ${PYTHON_GATO_FILES}

all: setup.py ghmm.py ${SWIG_INTERFACE_FILES}
	$(PYTHON) setup.py build

clean: setup.py
	$(PYTHON) setup.py clean
	rm -f message.txt

install: setup.py all
	$(PYTHON) setup.py install --prefix=${prefix}

check: ghmmunittests.py all
	PYTHONPATH=./build/lib.linux-i686-$(PYTHONVER)/:$(PYTHONPATH) LD_LIBRARY_PATH=../ghmm/.libs:$(LD_LIBRARY_PATH) $(PYTHON) ghmmunittests.py

coverage: ghmmunittests.py all
	coverage.py -e
	coverage.py -x ghmmunittests.py
	coverage.py -r ghmm.py

distclean: setup.py
	$(PYTHON) setup.py clean --all
	rm -f message.txt

.PHONY: coverage
