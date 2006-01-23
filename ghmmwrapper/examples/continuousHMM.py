import ghmm

# example code for a continuous HMM with gaussian emissions

F = ghmm.Float()  # emission domain of this model

A = [[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]]   # transition matrix
B = [[0.0,1.0],[-1.0,0.5], [1.0,0.2]]   # parameters of emission distributions in pairs of (mu, sigma)
pi = [1.0,0.0,0.0]   # initial probabilities per state

# generate model from parameters
model = ghmm.HMMFromMatrices(F,ghmm.GaussianDistribution(F), A, B, pi)

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
