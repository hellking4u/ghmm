#!/usr/bin/env python
#
# Package ghmm
#
#
#
from distutils.core import setup, Extension,DistutilsExecError
import os

# Adapted from Achim Gaedke's pygsl
def runtool(tool, argument_string):
    """ Run tool with the given arguments and return the first line (sans cr)
        of the output
    """
    command = os.popen(tool + ' ' + argument_string)
    output = command.readline()[:-1]
    command.close()
    if not output:
        raise DistutilsExecError, "could not run %s. Check your path." % tool
    return output

ghmmprefix  = runtool('ghmm-config' , '--prefix')
swiglib = runtool('swig','-swiglib')
swiglib_path  = os.path.split(swiglib)[0]
ghmmlib_path  = runtool('ghmm-config','--lib-prefix')
gslprefix = runtool('gsl-config', '--prefix')

print "********* PATHS ***********"
print "ghmmprefix " ,ghmmprefix
print "swiglib ", swiglib
print "swiglib_path ",swiglib_path
print "ghmmlib_path ", ghmmlib_path
print "**************************"


    
# BUG: Including 'ghmmwrapper.i' in Extension source list causes
# 'swig -python -o ghmmwrapper_wrap.c ghmmwrapper.i' to run.
# We just want: swig -c -python ghmmwrapper.i. Otherwise we get doubly
# defined symbols in link

print "================================================================================"
print "Please run the following command first: swig -noruntime -python -nodefault ghmmwrapper.i"
print "================================================================================"
   
setup(name="ghmmwrapper",
      version="0.6",
      description="Python Distribution Utilities",
      author="GHMM authors",
      author_email="ghmm@sf.net",
      url="http://ghmm.org",
      py_modules = ['ghmm','ghmmhelper','ghmmwrapper','modhmmer','xmlutil','DataStructures',
                    'Graph','GraphUtil', 'GatoGlobals','EditObjectAttributesDialog'],
      ext_modules = [Extension('_ghmmwrapper',
                               ['sdclass_change.c',
                                'gql.c', 'ghmmwrapper_wrap.c'],
                               include_dirs = [ghmmprefix + '/include', gslprefix + '/include'],
                               library_dirs = [ghmmlib_path ,swiglib_path],
                               libraries = ['gsl','stdc++','gsl','gslcblas','m','ghmm',
                               'swigpy' ] 
                               )
                     ]
     )

# EOF: setup.py
