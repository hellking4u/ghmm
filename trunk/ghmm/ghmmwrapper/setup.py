#!/usr/bin/env python
################################################################################
#
#       This file is part of the General Hidden Markov Model Library,
#       GHMM version __VERSION__, see http://ghmm.org
#
#       file:    setup.py
#       authors: Benjamin Georgi, Wasinee Rungsarityotin, Alexander Schliep
#
#       Copyright (C) 1998-2004 Alexander Schliep
#       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
#       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
#                               Berlin
#
#       Contact: schliep@ghmm.org
#
#       This library is free software; you can redistribute it and/or
#       modify it under the terms of the GNU Library General Public
#       License as published by the Free Software Foundation; either
#       version 2 of the License, or (at your option) any later version.
#
#       This library is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#       Library General Public License for more details.
#
#       You should have received a copy of the GNU Library General Public
#       License along with this library; if not, write to the Free
#       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#
#       This file is version $Revision$
#                       from $Date$
#             last change by $Author$.
#
################################################################################

from distutils.core import setup, Extension,DistutilsExecError
import os
import string

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
ghmmlib_libs  = string.split(runtool('ghmm-config','--libs'))[1:]

print "********* PATHS ***********"
print "ghmmprefix " ,ghmmprefix
print "swiglib ", swiglib
print "swiglib_path ",swiglib_path
print "ghmmlib_path ", ghmmlib_path
print "ghmmlib_libs ", ghmmlib_libs
print "**************************"


    
# BUG: Including 'ghmmwrapper.i' in Extension source list causes
# 'swig -python -o ghmmwrapper_wrap.c ghmmwrapper.i' to run.
# We just want: swig -c -python ghmmwrapper.i. Otherwise we get doubly
# defined symbols in link

print "================================================================================"
print "Please run the following command first: swig -python -nodefault ghmmwrapper.i"
print "================================================================================"
   
setup(name="ghmmwrapper",
      version="0.8",
      description="Python Distribution Utilities",
      author="GHMM authors",
      author_email="ghmm@sf.net",
      url="http://ghmm.org",
      py_modules = ['ghmm','ghmmhelper','ghmmwrapper','modhmmer','xmlutil','DataStructures',
                    'Graph','GraphUtil', 'GatoGlobals','EditObjectAttributesDialog','class_change'],
      ext_modules = [Extension('_ghmmwrapper',
                               ['sclass_change.c', 'pclasschange.c',
                                'gql.c', 'ghmmwrapper_wrap.c'],
                               include_dirs = [ghmmprefix + '/include'],
                               library_dirs = [ghmmlib_path, swiglib_path],
                               libraries = ['ghmm', 'm', 'pthread', 'xml2', 'z', 'swigpy'] #ghmmlib_libs.append('swigpy')
                               )
                     ]
     )

# EOF: setup.py
