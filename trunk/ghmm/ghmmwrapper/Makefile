#  author       : Wasinee Rungsarityotin and Benjamin Georgi
#  filename     : ghmmwrapper/Makefile
#  created      : DATE: September, 2003
#
#  $Id$
#
# __copyright__


#
# You must edit following variables:
#  TOP, SWIG to point to where SWIG is
#  PYTHON_INCLUDE, PYTHON_LIB
#

# Make sure these locate your installation of SWIG
TOP        = /project/algorithmics/
SWIG       = /project/algorithmics/bin/swig
TARGET     = ghmmwrapper
INTERFACE  = ghmmwrapper.i
SWIGOPT    = 


# These should locate your local installation of ghmm and ghmm++ (and xmlio if you use ghmm++)
GHMM_TOP   = /home/georgi/hmm/ghmm
GHMM_INCLUDE  = /home/georgi/hmm/ghmm
XMLIO_INCLUDE = /project/algorithmics/Sopra/xmlio-1.0.0/include
INCLUDE    =  -I$(GHMM_INCLUDE) -I$(XMLIO_INCLUDE)
LIBS       = -lstdc++ -L/usr/local/lib -lgsl -lgslcblas -lm -L$(GHMM_TOP)/lib -lghmm -lghmm++ -L/project/algorithmics/Sopra/xmlio-1.0.0/lib -lxmlio

#-lghmm++ -L/project/algorithmics/Sopra/xmlio-1.0.0/lib -lxmlio

#
# Make sure these locate your Python installation
#
PYTHON_INCLUDE= -I/package/python/2.3.0/linux/include/python2.3 -Ipackage/python/2.3.0/linux/lib/python2.3/config
PYTHON_LIB    = /package/python/2.3.0/linux/lib/python2.3/config

# Extra Python specific dynamic linking options
PYTHON_DLNK   =

#
# Must add the runtime library to LD_LIBRARY_PATH
#
all:
	$(SWIG) -c -python $(INTERFACE)
	$(CC) -c $(CFLAGS) sdclass_change.c $(INCLUDE)
	$(CC) -c $(CFLAGS) gql.c $(INCLUDE)
	$(CC) -c $(CFLAGS) read_cxml.c $(INCLUDE)
	$(CC) -c $(CFLAGS) ghmmwrapper_wrap.c $(INCLUDE) $(PYTHON_INCLUDE)
	$(CC) -shared sdclass_change.o gql.o read_cxml.o ghmmwrapper_wrap.o $(PYTHON_DLNK) $(LIBS) -o _ghmmwrapper.so -L$(TOP)/lib -lswigpy

clean:
	rm -f *_wrap* *.o *~ *.so mypython *.pyc .~*

