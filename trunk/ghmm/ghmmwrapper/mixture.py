#!/usr/bin/env python2.3
# 
# ghmm example: mixture.py
#
#
#
#
from ghmm import *
import numarray
import math

def sumlogs(a):
    """ Given a numarray.array a where a_i = log p_i, return log(sum p_i) """ 
    m = max(a)
    a = a - m
    result = 1.0
    for x in a:
        if x > -1.0e-16:
            result += math.exp(x) # XXX use approximation for speed
            
    result = math.log(result)
    result += m
    return result


hmmsFileName = 'test2.smo'
seqsFileName = 'test10.sqd'
maxiter = 100
eps = 0.1

models = HMMOpen.all(hmmsFileName)
print "# Read %d models from '%s'" % (len(models), hmmsFileName)
seqs = SequenceSet(Float(), seqsFileName)
print "# Read %d sequences from '%s'" % (len(seqs), seqsFileName)

done = 0
iter = 1
last_mixture_likelihood = -99999999.99
# The (nr of seqs x nr of models)-matrix holding the likelihoods
l = numarray.zeros((len(seqs), len(models)), numarray.Float)
logalpha = numarray.ones(len(models), numarray.Float) * math.log(1.0/len(models)) # Uniform alpha
print logalpha, numarray.exp(logalpha)
log_nrseqs = math.log(len(seqs))

while 1:
    # Score all sequences with all models
    for i, m in enumerate(models):
        loglikelihood = m.loglikelihoods(seqs) 
        # numarray slices: l[:,i] is the i-th column of l
        l[:,i] = numarray.array(loglikelihood)

    print l
    for i in xrange(len(seqs)):
        l[i] += logalpha # l[i] = ( log( a_k * P[seq i| model k]) )
    print l
    mixture_likelihood = numarray.sum(numarray.sum(l))
    print "# iter %s joint likelihood = %f" % (iter, mixture_likelihood) 

    improvement = mixture_likelihood - last_mixture_likelihood
    if iter > maxiter or improvement < eps:
        break

    # Compute P[model j| seq i]
    for i in xrange(len(seqs)):
        #l[i] += logalpha # l[i] = ( log( a_k * P[seq i| model k]) )
        # \sum_{k} a_k P[seq i| model k]
        seq_logprob = sumlogs(l[i])
        l[i] -= seq_logprob # l[i] = ( log P[model j | seq i] )
        
    # Compute priors alpha
    for i in xrange(len(models)):
        logalpha[i] = sumlogs(l[:,i]) - log_nrseqs

    print logalpha, numarray.exp(logalpha)
    
    for j, m in enumerate(models):
        # Set the sequence weight for sequence i under model m to P[m| i]
        for i in xrange(len(seqs)):
            seqs.setWeight(i,math.exp(l[i,j]))
        m.baumWelch(seqs, 10, 0.0001)

    iter += 1
