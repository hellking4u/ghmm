#/usr/bin/env python2.3
import ghmm

throws = IntegerRange(1,6)
dist = DiscreteDistiribution(throws)

edit_context = HMMEditingContext(throws, dist) # HMMEditContext
loaded = edit_context.AddState('loaded')
fair = edit_context.AddState('fair') # fair is a label used in pretty printing? XML only
fair.setEmissions([0.16, 0.16, 0.16, 0.16, 0.16, 0.16])
loaded.setEmissions([0.16, 0.9, 0.16, 0.16, 0.16, 0.16]) ## Auto-rescale?

edit_context.addTransition(fair, fair, 0.9) # or add_transition()
edit_context.addTransition(fair, loaded ,0.1) 
edit_context.addTransition(loaded, fair, 0.2)
edit_context.addTransition(loaded, loaded, 0.8)


lambda = edit_context() ## We cannot assure, when we delete states
# that initial indices remain valid. How can we find a state?


# How do we want to do that ...

observations = [...] # Or require EmissionSequence(throws,[...])
states = lamdba.viterbi(EmissionSequence(throws, observations)) 

# Alternatively: lambda knows EmissionDomain: Can call EmissionSequence
states = lamdba.viterbi(observations)

print #sequence + viterbi

# posterior decoding
forward = lamdba.forward(observations)

### How can we know, that the state with index 0 is 'fair' ???

posterior_fair = 


