from ghmm import *

print "*** EmissionSequence / SequenceSet ***"


HMMOpen = HMMOpenFactory(GHMM_FILETYPE_XML)
m4 = HMMOpen('/project/algorithmics/Sopra/python/simple2.xml')
print m4

alph = IntegerRange(0,7)

print "list input EmssionSequence:"
s = EmissionSequence(alph,[1,2,0,0,0,3,4])
print s



print "list(EmissionSequence)"
l = list(s)
print l



print "\nlist input SequenceSet:"
s2 = SequenceSet(alph,[ [1,2,3,4,5],[0,3,0],[4,3,2,2,1,1,1,1], [0,0,0,2,1],[1,1,1,1,1,1] ])
print s2

print "writing SequenceSet to file."
s2.write("blablatest.seq")

print "\nSequenceSet.getSubset"
s4 = s2.getSubset([0,2,4])
print s4


print "pointer input EmissionSequence"
s3 = s2[1] # call to SequenceSet.__getitem__
print s3


print "*** Discrete Emission Model ***"

m = HMMFromMatrices(DNA,DiscreteDistribution(DNA),
                       [[0.3,0.3,0.4],[0.6,0.1,0.3],[1.0,0.0,0.0]],
                       [[0.0,0.5,0.5,0.0],[0.1,0.0,0.8,0.1], [0.0,0.0,0.0,0.0]],
                       [1.0,0,0])

trans = m.getTransition(0,1)
print "a[0,1] = " + str(trans)

emission = m.getEmission(1)
print emission


#print m

#print "Sample:"					   
s4 = m.sample(4,15)
print str(s4) + "\n"

s5 = m.sampleSingle(10)
print str(s5) + "\n"

print "merging two sequences:"
s4.merge(s5)
print s4

print "training model"
m.baumWelch(s5)
#print m

print "Viterbi:"
path = m.viterbi(s5)
print str(path) + "\n"

print "forward"
logp1 = m.loglikelihood(s5)
print "logp = " + str(logp1) + "\n"

#print "logprob:"
#logp2 = m.logprob(s5,[0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1])
#print logp2, " -> " + str(2.71**logp2)

(alpha,scale) = m.forward(s5)
#print "alpha:\n" + str(alpha) + "\n"
#print "scale = " + str(scale) + "\n"	

beta = m.backward(s5,scale)
#print "beta = \n " + str(beta) + "\n"


print "\n\n\n *** Gaussian Model ***"
F = Float()            		   
m2 = HMMFromMatrices(F,GaussianDistribution(F),
                         [[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.0]],
                         [[0.0,1.0],[-1.0,0.5], [1.0,0.2]],
                         [1.0,0,0])

#print m2
trans = m2.getTransition(2,0)
print "a[2,2] = " + str(trans)
                         
#print "Sample:"
cs1 = m2.sample(4,15)                         
#print str(cs1) + "\n"

#print "SampleSingle:"
cs2 = m2.sampleSingle(10)                         
#print str(cs2) + "\n"

#print "Viterbi"
spath = m2.viterbi(cs1)
#print str(spath) + "\n"


#print "forward"
logp = m2.loglikelihood(cs1)    
#print "logp = " + str(logp) + "\n"

(salpha,sscale) = m2.forward(cs2)
#print "alpha:\n" + str(salpha) + "\n"
#print "scale = " + str(sscale) + "\n"	

beta = m2.backward(cs2,sscale)
#print "beta = \n " + str(beta) + "\n"


l = SequenceSetOpen(F,"seq_test.sqd")
#print l

print "*** Normalizing ***"
m3 = HMMFromMatrices(DNA,DiscreteDistribution(DNA),
                     [[10.0,10.0,10.0],[0.0,0.0,100.0],[25.0,25.0,50.0]],
                     [[10.0,0.0,010.0,0.0],[0.0,3.5,3.5,0.0], [5.0,5.0,5.0,5.0]],
                     [1.0,0,0])

#print m3

#m3.normalize()

#print m3

print "Writing to file:"
#m.write("er.log")
#m3.write("er.log")
#m2.write("er2.log")

#mList = [m,m2,m3]
#HMMwriteList("er.log",mList)


print m3
