from ghmm import *

sigma = IntegerRange(0,4)

m = HMMFromMatrices(sigma,DiscreteDistribution(sigma),
                       [[0.5,0.5,0.0,0.0],[0.0,0.0,1.0,0.0],[0.6,0.0,0.0,0.4],[1.0,0.0,0.0,0.0]],
                       [[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],
                       [1.0,0.0,0.0,0.0] )
#t2 = (0,-1,0,-1)
#m.setTieGroups(t2)


m2 = HMMFromMatrices(sigma,DiscreteDistribution(sigma),
                       [[0.3,0.3,0.4],[0.6,0.1,0.3],[0.3,0.3,0.4]],
                       [[0.3,0.2,0.3,0.2],[0.1,0.2,0.4,0.3], [0.25,0.25,0.25,0.25]],
                       [1.0,0,0])                       
                       
t = (0,-1,0)
m2.setTieGroups(t)

#m2.updateTieGroups()

print m2

train = m2.sample(50,1000)
#print m2
m2.baumWelch(train,3,10)
print m2
m2.updateTieGroups()
print "nach update"
print m2
m2.baumWelch(train,3,10)
