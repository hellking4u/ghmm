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

ghmmlib_path  = runtool('ghmm-config','--lib-prefix')

#ghmmprefix  = runtool('ghmm-config' , '--prefix')
#swiglib = runtool('swig','-swiglib')
#swiglib_path  = os.path.split(swiglib)[0]
#gslprefix = runtool('gsl-config', '--prefix')
#print "********* PATHS ***********"
#print "ghmmprefix " ,ghmmprefix
#print "swiglib ", swiglib
#print "swiglib_path ",swiglib_path
#print "ghmmlib_path ", ghmmlib_path
#print "**************************"

   
setup(name="Gato",
      version="x.xx",
      description="Python Distribution Utilities",
      author="Gato authors",
      py_modules = [ 'AnimatedAlgorithms', 'GatoFile', 'GraphEditor', 'ProbEditorDialogs', 
                     'AnimatedDataStructures', 'GatoGlobals', 'GraphUtil', 'ProbEditorWidgets',
                     'DataStructures', 'GatoIcons', 'Gred', 'TextTreeWidget',
                     'EditObjectAttributesDialog', 'GatoSystemConfiguration', 'HMMEd',
                     'TreeWidget', 'Embedder', 'GatoTest', 'HMMXML', 'convert-to-xml',
                     'GHMMXML', 'GatoUtil', 'MapEditor', 'logging', 'Gato', 'Graph',
                     'PlanarEmbedding', 'setup', 'Gato3D', 'GraphCreator', 'PlanarityTest',
                     'xmlutil', 'GatoConfiguration', 'GraphDisplay', 'ProbEditorBasics',
                     'GatoDialogs', 'GraphDisplay3D', 'ProbEditorContinuous'
                   ],
                   packages = ["editobj"]
     )

# EOF: setup




