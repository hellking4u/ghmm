from ghmm import *

sigma = IntegerRange(0,2)

m = HMMFromMatrices(sigma,DiscreteDistribution(sigma),
                       [[0.5,0.5],[1.0,0.0]],
                       [[1.0,0.0],[0.0,1.0]],
                       [1.0,0.0], labelList = ["GC rich","GC low"] )

print m.__class__.__name__
 
s = m.sample(3,40) 
print s

 
p = m.viterbiLabels(s)
print p
