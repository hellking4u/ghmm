from ghmm import *

print "*** EmissionSequence / SequenceSet ***"

alph = IntegerRange(0,7)

print "list input EmssionSequence:"
s = EmissionSequence(alph,[1,2,3,4])
print s

print "list input SequenceSet:"
s2 = SequenceSet(alph,[ [1,2,3,4,5],[0,3,0],[4,3,2,2,1,1,1,1] ])
print s2

print "pointer input EmissionSequence"
s3 = s2[1] # call to SequenceSet.__getitem__
print s3


print "*** Discrete Emission Model ***"

m = HMMFromMatrices(DNA,DiscreteDistribution(DNA),
                       [[0.0,1.0,0],[0.0,0.0,1.0],[1.0,0.0,0.0]],
                       [[1.0,0.0,0.0,0.0],[0.0,0.5,0.5,0.0], [0.0,0.0,0.0,0.0]],
                       [1.0,0,0])


#print m

print "Sample:"					   
s4 = m.sample(4,15)
print str(s4) + "\n"

s5 = m.sampleSingle(10)
print str(s5) + "\n"

print "Viterbi:"
path = m.viterbi(s5)
print str(path) + "\n"

print "forward"
logp1 = m.loglikelihood(s5)
print "logp = " + str(logp1) + "\n"

(alpha,scale) = m.forward(s5)
print "alpha:\n" + str(alpha) + "\n"
print "scale = " + str(scale) + "\n"	

beta = m.backward(s5,scale)
print "beta = \n " + str(beta) + "\n"


print "\n\n\n *** Gaussian Model ***"
F = Float()            		   
m2 = HMMFromMatrices(F,GaussianDistribution(F),
                         [[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.0]],
                         [[0.0,1.0],[-1.0,0.5], [1.0,0.2]],
                         [1.0,0,0])
                         
print "Sample:"
cs1 = m2.sample(4,15)                         
print str(cs1) + "\n"

print "SampleSingle:"
cs2 = m2.sampleSingle(10)                         
print str(cs2) + "\n"

print "Viterbi"
spath = m2.viterbi(cs1)
print str(spath) + "\n"


print "forward"
logp = m2.loglikelihood(cs1)    
print "logp = " + str(logp) + "\n"

(salpha,sscale) = m2.forward(cs2)
print "alpha:\n" + str(salpha) + "\n"
print "scale = " + str(sscale) + "\n"	

beta = m2.backward(cs2,sscale)
print "beta = \n " + str(beta) + "\n"
