################################################################################
#
#       This file is part of the General Hidden Markov Model Library,
#       GHMM version __VERSION__, see http://ghmm.org
#
#       file:    test_switch.py
#       authors: Benjamin Georgi
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
import ghmm
import ghmmwrapper
import ghmmhelper

A2 = [
      [[0.0,1.0,0.0],[0.0,0.0,1.0],[1.0,0.0,0.0]],
      [[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0]]
     ]
B2 = [[0.0,0.000001],[1.0,0.000001], [2.0,0.000001]]
pi2 = [1.0,0.0,0.0]
model2 = ghmm.HMMFromMatrices(ghmm.Float(),ghmm.GaussianDistribution(ghmm.Float), A2, B2, pi2)

ghmmwrapper.smodel_class_change_alloc(model2.cmodel)
ghmmwrapper.setPythonSwitching(model2.cmodel,"class_change","getClass")
#ghmmwrapper.setSwitchingFunction(model2.cmodel)

#print model2
seq = model2.sample(2,30)
