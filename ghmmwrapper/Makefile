TOP        = /project/algorithmics/Sopra/SWIG-1.3.17
SWIG       = $(TOP)/swig
TARGET     = ghmmwrapper
INTERFACE  = ghmmwrapper.i
SWIGOPT    = 
# GHMM_TOP   = /project/algorithmics/Sopra/ghmm-0.3.0
# GHMM_TOP   = /project/algorithmics/Sopra/ghmm
GHMM_TOP = /home/georgi/bachelor/ghmm
GHMM_INCLUDE  = /home/georgi/bachelor/ghmm
XMLIO_INCLUDE = /project/algorithmics/Sopra/xmlio-1.0.0/include
INCLUDE    =  -I$(GHMM_INCLUDE) -I$(XMLIO_INCLUDE)
LIBS       = -lstdc++ -L/usr/local/lib -lgsl -lgslcblas -lm -L$(GHMM_TOP)/lib -lghmm -lghmm++ -L/project/algorithmics/Sopra/xmlio-1.0.0/lib -lxmlio
CFLAGS = -g

# Make sure these locate your Python installation
# PYTHON_INCLUDE= -DHAVE_CONFIG_H -I/home/rungsari/Python-2.2/include/python2.2 -I/home/rungsari/Python-2.2/lib/python2.2/config
PYTHON_INCLUDE= -I/package/python/2.3.0/Linux/include/python2.3/ -Ipackage/python/2.3.0/Linux/include/python2.3/config
PYTHON_LIB    = /package/python/2.3.0/Linux/lib/python2.3/config/ 


# Extra Python specific dynamic linking options
PYTHON_DLNK   =

#
# Must add the runtime library to LD_LIBRARY_PATH
#
all:
	$(SWIG) -c -python $(INTERFACE)
	$(CC) -c $(CFLAGS) read_cxml.c $(INCLUDE)
	$(CC) -c $(CFLAGS) gql.c $(INCLUDE)
	$(CC) -c $(CFLAGS) ghmmwrapper_wrap.c $(INCLUDE) $(PYTHON_INCLUDE)
	$(CC) -shared gql.o read_cxml.o ghmmwrapper_wrap.o $(PYTHON_DLNK) $(LIBS) -o _ghmmwrapper.so -L$(TOP)/lib -lswigpy

clean:
	rm -f *_wrap* *.o *~ *.so mypython *.pyc .~*
