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




 # Interpretation of B matrix for the mixture case (Example with three states and two components each):
                #  B = [ 
                #      [ ["mu11","mu12"],["sig11","sig12"],["w11","w12"]   ],
                #      [  ["mu21","mu22"],["sig21","sig22"],["w21","w22"]  ],
                #      [  ["mu31","mu32"],["sig31","sig32"],["w31","w32"]  ],
                #      ]



F = ghmm.Float()            		   
m2 = ghmm.HMMFromMatrices(F,ghmm.GaussianMixtureDistribution(F),
                         [[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]],
                         
                         [
                         [ [0.0,1.0],[1.0,0.5], [0.6,0.4] ],
                         [ [2.0,1.0],[1.0,0.5], [0.8,0.2] ],
                         [ [15.0,1.0],[1.0,0.5], [0.8,0.2] ]
                         ],
                         
                         [1.0,0,0])





A = [[0.9,0.1],[0.2,0.8]]
B = [ [[0.0,0.0],[0.1,0.1],[0.5,0.5]],
    [[10.0,10.0],[0.1,0.1],[0.5,0.5]] ]
pi = [0.0,1.0]
model = ghmm.HMMFromMatrices(ghmm.Float(),ghmm.GaussianMixtureDistribution(ghmm.Float), A, B, pi)

#short = ghmm.EmissionSequence(ghmm.Float(),[0.0,10.0,10.0,0.0,10.0,0.0,10.0,10.0,0.0])
#short2 = ghmm.EmissionSequence(ghmm.Float(),[10.0,0.0]*6)


#print model
#model.baumWelch(short2,1,0.2)
#print model

#def getClass(seq, k,t):
#    classes = [0] * 10 + [1] * 10
#    return classes[t]


swA = [
    [[0.9,0.1],[0.2,0.8]],
    [[0.3,0.7],[0.7,0.3]]
    ]

swB = [ [[0.0,0.0],[0.1,0.1],[0.5,0.5]],
    [[10.0,10.0],[0.1,0.1],[0.5,0.5]] ]
swpi = [0.5,0.5]
swmodel = ghmm.HMMFromMatrices(ghmm.Float(),ghmm.GaussianMixtureDistribution(ghmm.Float), swA, swB, swpi)

ghmmwrapper.smodel_class_change_alloc(swmodel.cmodel)
ghmmwrapper.setPythonSwitching(swmodel.cmodel,"class_change","getClass")   # "class_change"

swshort = ghmm.EmissionSequence(ghmm.Float(),[10.0]+[0.0]*5+[0.0]+[10.0]*5)

print swshort

print swmodel
swmodel.baumWelch(swshort,5,0.2)
print swmodel













#A = [ [[0.0,1.0],[0.0,0.0]],
#      [[1.0,0.0],[1.0,0.0]]
#    ]


                
#A2 = [ [[0.0,1.0,0.0],[0.0,0.0,1.0],[1.0,0.0,0.0]],
#      [[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0]] ]
     
#B2 = [ [[0.0,0.000001],[0.4,0.4],[0.5,0.5]],
#    [[1.0,0.000001],[0.3,0.3],[0.5,0.5]],
#     [[2.0,0.000001],[0.8,0.8],[0.5,0.5]] ]


#pi2 = [1.0,0.0,0.0]
#model2 = ghmm.HMMFromMatrices(ghmm.Float(),ghmm.GaussianMixtureDistribution(ghmm.Float), A2, B2, pi2)

#ghmmwrapper.smodel_class_change_alloc(model2.cmodel)
#ghmmwrapper.setPythonSwitching(hmm.cmodel,"class_change","getClass")

#ghmmwrapper.setPythonSwitching(model2.cmodel,"class_change","getClass")
#ghmmwrapper.setSwitchingFunction(model2.cmodel)

#print model2
#seq = model2.sample(2,30)

#seq = EmissonSequence(ghmm.Float(),



#model2.write("test_switch.smo")
#m2 = ghmm.HMMOpen("test_switch.smo")
#print m2
