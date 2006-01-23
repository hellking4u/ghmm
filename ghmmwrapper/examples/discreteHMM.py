import ghmm

# example code for a discrete HMM for generating random DNA sequences

A = [[0.3,0.3,0.4],[0.6,0.1,0.3],[1.0,0.0,0.0]]  # transition matrix
B = [[0.2,0.3,0.3,0.2],[0.1,0.2,0.6,0.1], [0.25,0.25,0.25,0.25]]  # emission matrix
pi = [1.0,0.0,0.0]  # initial probabilities

# generate model from parameter matrices
model = ghmm.HMMFromMatrices(ghmm.DNA,ghmm.DiscreteDistribution(ghmm.DNA), A, B, pi)

# sample single sequence of length 50
seq = model.sampleSingle(50)

# sample 10 sequences of length 50
seq_set = model.sample(10,50)

# get log P(seq | model)
logp = model.loglikelihood(seq)

# cacluate viterbi path 
path = model.viterbi(seq)

# train model parameters
model.baumWelch(seq_set,5,0.01)
