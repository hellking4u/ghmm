#!/usr/bin/env python
#
# Package pyghmm
#
#
#
from distutils.core import setup, Extension
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

xmlioprefix = runtool('xmlio-config', '--prefix')
ghmmprefix  = runtool('ghmm-config' , '--prefix')
swiglib = runtool('swig','-swiglib')
swiglib_path  = os.path.split(swiglib)[0]
xmliolib_path = xmlioprefix + '/lib'
ghmmlib_path  = runtool('ghmm-config','--lib-prefix')
    
# BUG: Including 'ghmmwrapper.i' in Extension source list causes
# 'swig -python -o ghmmwrapper_wrap.c ghmmwrapper.i' to run.
# We just want: swig -c -python ghmmwrapper.i. Otherwise we get doubly
# defined symbols in link

print "================================================================================"
print "BUG:  Please run the following command: swig -c -python ghmmwrapper.i"
print "================================================================================"
   
setup(name="ghmmwrapper",
      version="0.6",
      description="Python Distribution Utilities",
      author="GHMM authors",
      author_email="ghmm@sf.net",
      url="http://ghmm.org",
      py_modules = ['ghmm','ghmmwrapper','modhmmer'],
      ext_modules = [Extension('_ghmmwrapper',
                               ['sdclass_change.c',
                                'gql.c', 'read_cxml.c', 'ghmmwrapper_wrap.c'],
                               include_dirs = [xmlioprefix + '/include',ghmmprefix + '/include'],
                               library_dirs = [xmliolib_path , ghmmlib_path ,swiglib_path],
                               libraries = ['gsl','stdc++','gsl','gslcblas','m','ghmm',
                                            'ghmm++', 'xmlio', 'swigpy'
                                            ]
                               )
                     ]
     )

# EOF: setup.py
