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

from distutils.core import setup,Extension

# BUG: Including 'ghmmwrapper.i' in Extension source list causes
# 'swig -python -o ghmmwrapper_wrap.c ghmmwrapper.i' to run.
# We just want: swig -c -python ghmmwrapper.i. Otherwise we get doubly
# defined symbols in link
# NOT A BUG ANYMORE? I get only doubly defined symbols if the swig interface file
#                    and the generated c file are in the extension list.
#                    I think that is expected behaviour.
#                    Tested with swig 1.3.19, 1.3.21 and 1.3.25 
   
setup(name="ghmmwrapper",
      version="0.8",
      description="Python Distribution Utilities",
      author="GHMM authors",
      author_email="ghmm-list@lists.sourceforge.net",
      url="http://ghmm.org",
      py_modules  = ['ghmm','ghmmhelper','ghmmwrapper','modhmmer','class_change'],
      packages    = ['ghmm_gato'],
      ext_modules = [Extension('_ghmmwrapper',
                               ['sclass_change.c', 'pclasschange.c', 'gql.c', 'ghmmwrapper.i'],
                               include_dirs = ['..'],
                               library_dirs = ['../ghmm/.libs'],
                               libraries = ['ghmm', 'm', 'pthread', 'xml2', 'z'],
                               extra_compile_args = ["-O2", "-pipe", "-Wall"], # -g might help debugging
                               depends = ['wrapper_alphabet.i', 'wrapper_cmodel.i', 'wrapper_cseq.i',
                                          'wrapper_dmodel.i', 'wrapper_dpmodel.i', 'wrapper_dpseq.i',
                                          'wrapper_dseq.i', 'wrapper_xmlfile.i']
                               )
                     ]
     )

# EOF: setup.py
