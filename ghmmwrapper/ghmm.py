#!/usr/bin/env python2.3
################################################################################
#
#       This file is part of GHMM (General Hidden Markov Model library) 
#
#       file:    ghmm.py
#       authors: Benjamin Georgi, Wasinee Rungsarityotin, Alexander Schliep
#
#       Copyright (C) 2003-2004, Alexander Schliep and MPI Molekulare Genetik, Berlin
#                                   
#       Contact: schliep@molgen.mpg.de         
#
#       Information: http://ghmm.org
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
#       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
#
#
#       This file is version $Revision$ 
#                       from $Date$
#             last change by $Author$.
#
################################################################################

"""

The Design of ghmm.py 

HMMs are stochastic models which encode a probability density over
sequences of symbols. These symbols can be discrete letters (A,C,G and
T for DNA; 1,2,3,4,5,6 for dice), real numbers (weather measurement
over time: temperature) or vectors of either or the combination
thereof (weather again: temperature, pressure, percipitation).

Note: We will always talk about emissions, emission sequence and so
forth when we refer to the sequence of symbols. Another name
for the same object is observation resp. observation sequence.


The objects one has to deal with in HMM modelling are the following

1) The domain the emissions come from: the EmissionDomain. Domain
   is to be understood mathematically and to encompass both discrete,
   finite alphabets and fields such as the real numbers or intervals
   of the reals.

   For technical reasons there can be two representations of an emission
   symbol: an external and an internal. The external representation
   is the view of the application using ghmm.py. The internal one
   is what is used in both ghmm.py and the ghmm C-library. Representations
   can coincide, but this is not guaranteed. Most prominently discrete
   alphabets of size k are represented as [0,1,2,...,k] internally.
   It is the domain objects job to provide a mapping between representations
   in both directions.

   NOTE: Do not make assumptions about the internal representations. It might
   change.
   

2) Every domain has to afford a distribution, which is usually
   parameterized. A distribution associated with a domain
   should allow us to compute $\Prob[x| distribution parameters]$
   efficiently.

   The distribution defines the *type* of distribution which
   we will use to model emissions in *every state* of the HMM.
   The *type* of distribution will be identical for all states,
   their *parameterizations* will differ from state to state. 
   

3) We will consider a Sequence of emissions from the same emission
   domain and very often sets of such sequences: SequenceSet


4)
   The HMM: The HMM consists of two major components: A Markov chain
   over states (implemented as a weighted directed graph with
   adjacency and inverse-adjacency lists) and the emission
   distributions per-state. For reasons of efficiency the HMM itself
   is *static*, as far as the topology of the underlying Markov chain
   (and obviously the EmissionDomain) are concerned. You cannot add or
   delete transitions in an HMM.

   Transition probabilities and the parameters of the per-state
   emission distributions can be easily modified. Particularly,
   Baum-Welch reestimation is supported.  While a transition cannot be
   deleted from the graph, you can set the transition probability to
   zero, which has the same effect from the theoretical point of
   view. However, the corresponding edge in the graph is still
   traversed in the computation. 

   States in HMMs are referred to by their integer index. State sequences
   are simply list of integers.

   If you want to store application specific data for each state you have to
   do it yourself.

   Subclasses of HMM implement specific types of HMM. The type depends
   on the EmissionDomain, the Distribution used, the specific
   extensions to the 'standard' HMMs and so forth

   
5) HMMFactory: This provides a way of constucting HMMs. Classes derived
   from HMMFactory allow to read HMMs from files, construct them explicitly
   from, for a discrete alphabet, transition matrix, emission matrix and prior
   or serve as the basis for GUI-based model building.

   There are two ways of using the HMMFactory.

   Static construction:

   HMMOpen(fileName) # Calls an object of type HMMOpen instantiated in ghmm
   HMMOpen(fileName, type=HMM.FILE_XML)
   HMMFromMatrix(emission_domain, distribution, A, B, pi) # B is a list of distribution parameters

   # XXX do we need dsitribution here?

   Dynamic construction: Providing a context for dynamically
   editing existing HMMs or creating them from scratch

   HMMEditingContext() # Create an object
   HMMEditingContext(hmm) # Create an object from existing HMM
   
   Examples:

   hmm = HMMOpen('some-hmm.xml')

   hmm_context = HMMEditingContext(hmm) # just reads hmm
   hmm_context.addTransition(4,5, 0.3) # normalization will occurr
   hmm_context.addTransition(5,6, 0.1)

   hmm = hmm_context() # Creates a new hmm 

   hmm.bla ....

   

"""


import ghmmwrapper
import ghmmhelper

#-------------------------------------------------------------------------------
#- Exceptions ------------------------------------------------------------------

class GHMMError(Exception):
    """Base class for exceptions in this module."""

#    def __init__(self, expression, message):
#        self.expression = expression
#        self.message = message

class UnknownInputType(GHMMError):
    def __init__(self):
        pass

class NoValidCDataType(GHMMError):
    def __init__(self):   
        pass
#-------------------------------------------------------------------------------
#- EmissionDomain and derived  -------------------------------------------------
class EmissionDomain:
    """ Abstract base class for emissions produced by an HMM.

        There can be two representations for emissions:
        1) An internal, used in ghmm.py and the ghmm C-library
        2) An external, used in your particular application

        Example:
        The underlying library represents symbols from a finite,
        discrete domain as integers (see Alphabet).

        EmissionDomain is the identity mapping
    """

    def internal(self, emission):
        """ Given a emission return the internal representation
        """
        return emission


    def internalSequence(self, emissionSequence):
        """ Given a emissionSequence return the internal representation
        """
        return emissionSequence


    def external(self, internal):
        """ Given an internal representation return the
            external representation
        """
        return internal

    def externalSequence(self, internalSequence):
        """ Given a sequence with the internal representation return the
            external representation
        """
        return internalSequence


    def isAdmissable(self, emission):
        """ Check whether emission is admissable (contained in) the domain
            raises GHMMOutOfDomain else
        """
        return None


class Alphabet(EmissionDomain):
    """ Discrete, finite alphabet
    
    """
    def __init__(self, listOfCharacters):
        """ Creates an alphabet out of a listOfCharacters """
        self.listOfCharacters = listOfCharacters
        self.index = {} # Which index belongs to which character
        i = 0
        for c in self.listOfCharacters:
            self.index[c] = i
            i += 1
        self.CDataType = "int" # flag indicating which C data type should be used

    def internal(self, emission):
        """ Given a emission return the internal representation
        """
        return self.index[emission]


    def internalSequence(self, emissionSequence):
        """ Given a emission_sequence return the internal representation
        """
        result = copy.deepcopy(emissionSequence)
        for i in xrange(len(result)):
            result[i] = self.index[result[i]]
        return result


    def external(self, internal):
        """ Given an internal representation return the
            external representation
        """
        return self.listOfCharacters[internal]

    def externalSequence(self, internalSequence):
        """ Given a sequence with the internal representation return the
            external representation
        """
        result = copy.deepcopy(internalSequence)
        for i in xrange(len(result)):
            result[i] = self.listOfCharacters[result[i]]
        return result

    def isAdmissable(self, emission):
        """ Check whether emission is admissable (contained in) the domain
            raises GHMMOutOfDomain else
        """
        return emission in self.listOfCharacters

    def size(self):
        return len(self.listOfCharacters)


DNA = Alphabet(['A','C','G','T'])
AminoAcids = Alphabet(['A','C','D','E','F','G','H','I','K','L',
                       'M','N','P','Q','R','S','T','V','W','Y'])
def IntegerRange(a,b):
    return Alphabet(range(a,b))


class Float(EmissionDomain):

    def __init__(self):
        self.CDataType = "double" # flag indicating which C data type should be used

    def isAdmissable(self, emission):
        """ Check whether emission is admissable (contained in) the domain
            raises GHMMOutOfDomain else
        """
        return isinstance(a,float)
   

    
#-------------------------------------------------------------------------------
#- Distribution and derived  ---------------------------------------------------
class Distribution:
    """ Abstract base class for distribution over EmissionDomains
    """

    # From Spass (S. Rahmann):
    # add density, mass, cumuliative dist, quantils, sample, fit pars,
    # moments


class DiscreteDistribution(Distribution):
    """ A DiscreteDistribution over an Alphabet: The discrete distribution
        is parameterized by the vectors of probabilities.


    """

    def __init__(self, alphabet):
        self.alphabet = alphabet

        

class ContinousDistribution(Distribution):
    pass

class GaussianDistribution(ContinousDistribution):

    def __init__(self, domain):
        self.emissionDomain = domain
        #self.mu = mu
        #self.sigma2 = sigma2

    def set(self, (mu, sigma2)):
        self.mu = mu
        self.sigma2 = sigma2

    def get(self):
        return (self.mu, self.sigma2)
        

class MixtureContinousDistribution(Distribution):
    pass


#-------------------------------------------------------------------------------
#- Sequence, SequenceSet and derived  ------------------------------------------
class EmissionSequence(list):
    """ An EmissionSequence contains the *internal* representation of
        a sequence of emissions. It also contains a reference to the
        domain where the emission orginated from.
    """

    def __init__(self, emissionDomain, sequenceInput ,inputType = "pylist"):
        self.emissionDomain = emissionDomain # XXX identical name problematic ?
        
        if self.emissionDomain.CDataType == "int": # underlying C data type is integer
            if inputType == "pylist":
                l = len(sequenceInput)
                seq = ghmmhelper.list2arrayint(sequenceInput)
                self.cseq = sequence_t()
                self.cseq.seq
                
            elif inputType == "fromFile":    
                pass
            elif inputType == "fromCpointer": # for internal use mostly
                pass
            else:    
                raise UnknownInputType, "inputType " + str(inputType) + " not recognized."
        
        elif self.emissionDomain.CDataType == "double": # underlying C data type is double
            if inputType == "pylist":
                pass
            elif inputType == "fromFile":    
                pass
            elif inputType == "fromCpointer": # for internal use mostly
                pass
            else:    
                raise UnknownInputType, "inputType " + str(inputType) + " not recognized."
        
        
        else:
            raise NoValidCDataType, "C data type " + str(self.emissionDomain.CDataType) + " invalid."
        
        
        self.c_seq 
        
#    def sequenceSet(self):
#        """ Make a one-element SequenceSet out of me """

class SequenceSet:
    def __init__(self, emissionDomain, fromList = None):
        self.emissionDomain = emissionDomain # XXX identical name problematic ?
        self.c_seq 
    

def SequenceSetOpen(fileName):
    c_ptr = ghmmwrapper.seq_d_read(fileName)  
    return SequenceSet(c_ptr)



#-------------------------------------------------------------------------------
#- HMMFactory and derived  -----------------------------------------------------
class HMMFactory:
    """ A HMMFactory is the base class of HMM factories.
        A HMMFactory has just a constructor and a () method
    """


GHMM_FILETYPE_SMO = 'smo'
GHMM_FILETYPE_XML = 'xml'

class HMMOpenFactory(HMMFactory):

    def __init__(self, defaultFileType=None):
        if defaultFileType:
            self.defaultFileType = defaultFileType

    def __call__(self, fileName, fileType=None):
        hmmClass = self.determineHMMClass("bla.xml") # Return proper class
        #m = hmmClass("bla.xml")
        #return m

    def determineHMMClass(self, fileName):
        pass 

HMMOpen = HMMOpenFactory(GHMM_FILETYPE_SMO)

class HMMFromMatricesFactory(HMMFactory):
    def __call__(self, emissionDomain, distribution, A, B, pi):

        if isinstance(emissionDomain,Alphabet):

            if isinstance(distribution,DiscreteDistribution):
                # HMM has discrete emissions over finite alphabet: DiscreteEmissionHMM
                cmodel = ghmmwrapper.model()

                cmodel.N = len(A)
                cmodel.M = emissionDomain.size()
                cmodel.prior = -1 # No 
                cmodel.name = 'Unused'
                
                states = ghmmwrapper.arraystate(cmodel.N)
                
                #initialize states
                for i in range(cmodel.N):
                    state = ghmmwrapper.get_stateptr(states,i)
                    state.b = ghmmhelper.list2arrayd(B[i])
                    state.pi = pi[i]
                    
                    #if (sum(B[i]) <= epsilon):
                    #        silent_states.append(1)
                    #        silent_flag = 4
                    #else:
                    #        silent_states.append(0)

                    #set out probabilities
                    
                    trans = ghmmhelper.extract_out(A[i])
                    state.out_states = trans[0]
                    state.out_id = trans[1]
                    state.out_a = trans[2]
                    #set "in" probabilities
                    # XXX Check whether A is numarray or Python
                    A_col_i = map( lambda x: x[i], A)
                    # Numarray use A[,:i]
                    trans = ghmmhelper.extract_out(A_col_i)
                    state.in_states = trans[0]
                    state.in_id = trans[1]
                    state.in_a = trans[2]
                    #fix probabilities by reestimation, else 0
                    state.fix = 0
                    
                cmodel.s = states
                #model.silent = int_array(model.N)
                #model.model_type = silent_flag
                #plist2intarray(model.silent, silent_states, model.N)
                return DiscreteEmissionHMM(emissionDomain, distribution, cmodel)

            else:
                raise GHMMError(type(distribution), "Not a valid distribution for Alphabet") 
        else:
			
            if isinstance(distribution,GaussianDistribution):
				
				cmodel = ghmmwrapper.smodel()
				cmodel.N = len(A)
				cmodel.M = 1 # Number of mixture componenent for emission distribution
				cmodel.prior = -1 # Unused
				cmodel.cos = 1  # number of transition classes in GHMM
				states = ghmmwrapper.arraysstate(cmodel.N)

				# XXX ? switching function ? XXX
				# Switching functions and transition classes are handled
				# elswhere

				#initialize states
				for i in range(cmodel.N):
					state = ghmmwrapper.get_sstate_ptr(states,i)
					state.pi = pi[i]

					# allocate arrays of emmission parameters
					state.c = ghmmhelper.list2arrayd([1.0]) # Mixture weights. Unused
					(mu, sigma) = B[i]
					state.mue = ghmmhelper.list2arrayd([mu]) #mu = mue in GHMM C-lib.
					state.u = ghmmhelper.list2arrayd([sigma])

					#set out probabilities
					trans = ghmmhelper.extract_out_probs([A[i]], cmodel.cos) # cos = 1
					state.out_states = trans[0]
					state.out_id = trans[1]
					state.out_a = trans[2].array 

					#set "in" probabilities
					A_col_i = map( lambda x: x[i], A)
					trans = ghmmhelper.extract_out_probs([A_col_i],cmodel.cos) # cos = 1
					state.in_states = trans[0]
					state.in_id = trans[1]
					state.in_a = trans[2].array
				
					state.fix = 0 # if fix = 1, exclude state's probabilities from reestimation

				#append states to model
				cmodel.s = states
				return GaussianEmissionHMM(emissionDomain, distribution, cmodel)

            else:
                raise GHMMError(type(distribution),
								"Cannot construct model for this domain/distribution combination") 


HMMFromMatrices = HMMFromMatricesFactory()



#-------------------------------------------------------------------------------
#- HMM and derived  ------------------------------------------------------------
class HMM:

    def __init__(self, emissionDomain, distribution, cmodel):
        self.emissionDomain = emissionDomain
        self.distribution = distribution
        self.cmodel = cmodel


    def loglikelihood(self, emissionSequences):
        """ Compute log( P[emissionSequences| model]) using the forward algorithm

            emission_sequences can either be a SequenceSet or a Sequence

            Result: log( P[emissionSequences| model]) of type float
                    numarray vector of floats
            
            Note: The implementation will not compute the full forward matrix 
        """
        log_p = double_array(1)		
        i = self.cmodel.N
        logp_sum = 0
        for seq_nr in range(mysequence.seq_c.seq_number):
            t = get_arrayint(mysequence.seq_c.seq_len,seq_nr)
            alpha = matrix_d_alloc(t,i)
            scale = double_array(t)
            
            seq = mysequence[seq_nr]
            error = func(self.model, seq.seq,seq.length, alpha, scale, log_p)
            if error == -1:
                print "ERROR: Forward finished with -1: Sequence " + str(seq_nr) + " cannot be build."
                exit
            logp_sum  += get_arrayd(log_p,0)
            print "Seq " + str(seq_nr) + ": " + str(get_arrayd(log_p,0))
            
        return (logp_sum,scale) # XXX TEST



    ## Further Marginals ...



    def logprob(self, emissionSequence, stateSequence):
        """ log P[ emissionSequence, stateSequence| m] """
        pass
        

    def baumWelch(self, trainingSequences, nrSteps, loglikelihoodCutoff):
        """ Reestimate the model parameters given the training_sequences.
            Perform at most nr_steps until the improvement in likelihood
            is below likelihood_cutoff
        
            training_sequences can either be a SequenceSet or a Sequence
  
            Result: Final loglikelihood
        """
        self.baumWelchSetup(trainingSequences)
        (steps_made, loglikelihood_array, scale_array) = self.baumWelchStep(nrSteps,
                                                                            loglikelihoodCutoff)
        return loglikelihood_array[-1]

    def baumWelchSetup(self, trainingSequences):
        """ Setup necessary temporary variables for Baum-Welch-reestimation.
            Use baum_welch_setup and baum_welch_step if you want more control
            over the training, compute diagnostics or do noise-insertion

            training_sequences can either be a SequenceSet or a Sequence
        """
        pass
    
    def baumWelchStep(self, nrSteps, loglikelihoodCutoff):
        """ Setup necessary temporary variables for Baum-Welch-reestimation.
            Use baum_welch_setup and baum_welch_step if you want more control
            over the training, compute diagnostics or do noise-insertion

            training_sequences can either be a SequenceSet or a Sequence
        """
        pass
    
    def baumWelchDelete(self):
        """ Delete the necessary temporary variables for Baum-Welch-reestimation """
        # Needed ?
        pass

    def forward(self, emissionSequence):
        """

            Result: the (N x T)-matrix containing the forward-variables
                    and the scaling vector
        """
        pass

    def backward(self, emissionSequence, scalingVector):
        """

            Result: the (N x T)-matrix containing the backward-variables
        """

    

    def viterbi(self, emissionSequences):
        """ Compute the Viterbi-path for each sequence in emission_sequences

            emission_sequences can either be a SequenceSet or a Sequence

            Result: [q_0, ..., q_T] the viterbi-path if emission_sequences is a Sequence
                    [[q_0^0, ..., q_T^0], ..., [q_0^k, ..., q_T^k]} for a k-sequence
                    SequenceSet
        """

    
    def sample(self, T):
        """ Sample emission sequences 


        """

    def sampleSingle(self, T):
        """ Sample a single emission sequence of length at most T.
            Returns a Sequence object.
        """

    def state(self, stateLabel):
        """ Given a stateLabel return the integer index to the state """
        pass

    def getInitial(self, i):
        """ Accessor function for the initial probability \pi_i """
        pass

    def setInitial(self, i, j, prob):
        """ Accessor function for the initial probability \pi_i """
        pass

    def getTransition(self, i, j):
        """ Accessor function for the transition a_ij """
        pass

    def setTransition(self, i, j, prob):
        """ Accessor function for the transition a_ij """
        pass

    def getEmission(self, i):
        """ Accessor function for the  """
        pass

    def setEmission(self, i, distributionParemters):
        """ Set the emission distribution parameters """
        pass

    def normalize(self):
        """ Normalize transition probs, emission probs (if applicable) """
        pass
        
    def randomize(self, noise_level):
        """ """
        pass
    

class DiscreteEmissionHMM(HMM):
    
    def __init__(self, emissionDomain, distribution, cmodel):
        HMM.__init__(self, emissionDomain, distribution, cmodel)


    def __str__(self):
        hmm = self.cmodel
        strout = "\nOverview of HMM:\n"
        strout += "\nNumber of states: "+ str(hmm.N)
        strout += "\nSize of Alphabet: "+ str(hmm.M)

        for k in range(hmm.N):
            state = ghmmwrapper.get_stateptr(hmm.s, k)
            strout += "\n\nState number "+ str(k) +":"
            strout += "\nInitial probability: " + str(state.pi)
            #strout += "\nsilent state: " + str(get_arrayint(self.model.silent,k))
            strout += "\nOutput probabilites: "
            for outp in range(hmm.M):
                strout+=str(ghmmwrapper.get_arrayd(state.b,outp))+", "
            strout += "\nOutgoing transitions:"
            for i in range( state.out_states):
                strout += "\ntransition to node " + str( ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.out_a,i))
            strout +=  "\nIngoing transitions:"
            for i in range(state.in_states):
                strout +=  "\ntransition from node " + str( ghmmwrapper.get_arrayint(state.in_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.in_a,i))
                strout += "\nint fix:" + str(state.fix) + "\n"
            #strout += "Silent states: \n"
            #for k in range(hmm.N):
            #strout += str(get_arrayint(self.model.silent,k)) + ", "
            #strout += "\n"
        return strout
    


class GaussianEmissionHMM(HMM):
    """ GaussianEmissionHMM are HMMs which have a Gaussian distribution
        per state determining the emission probababilities
    """
  
    def __init__(self, emissionDomain, distribution, cmodel):
        HMM.__init__(self, emissionDomain, distribution, cmodel)
    

    def loglikelihood(self, emissionSequences):
        """ Compute log( P[emissionSequences| model]) using the forward algorithm

            emission_sequences can either be a SequenceSet or a Sequence

            Result: log( P[emissionSequences| model]) of type float
                    numarray vector of floats
            
            Note: The implementation will not compute the full forward matrix 
        """

    def loglikelihood_sqd(self, sequenceSet):
        # XXX REMOVE soon XXX
        """ Compute log( P[emissionSequences| model]) using the forward algorithm

            emission_sequences can either be a SequenceSet or a Sequence

            Result: log( P[emissionSequences| model]) of type float
                    numarray vector of floats
            
            Note: The implementation will not compute the full forward matrix 
        """
        likelihood = ghmmwrapper.double_array(sequenceSet.seq_number)
        result = ghmmwrapper.smodel_individual_likelihoods(self.cmodel, sequenceSet, likelihood)
        return likelihood # Caller owns it 
	
    def getTransition(self, i, j):
        """ Accessor function for the transition a_ij """
        return ghmmwrapper.smodel_get_transition(self.cmodel, i, j, 0)

    def setTransition(self, i, j, prob):
        """ Accessor function for the transition a_ij """
        ghmmwrapper.smodel_set_transition(self.cmodel, i, j, 0, float(prob))
       
    def getEmission(self, i):
        """ Return (mu, sigma^2)  """
        state = ghmmwrapper.get_sstate(self.cmodel, i)
        mu = ghmmwrapper.get_arrayd(state.mue, 0)
        sigma = ghmmwrapper.get_arrayd(state.u,0)
        return (mu, sigma)
        
    def setEmission(self, i, (mu, sigma)):
        """ Set the emission distributionParameters for state i """
        state = ghmmwrapper.get_sstate(self.cmodel, i)
        ghmmwrapper.set_arrayd(state.mue, 0, float(mu))  # GHMM C is german: mue instead of mu 
        sigma = ghmmwrapper.set_arrayd(state.u, 0, float(sigma))       
   

    def __str__(self):
        # XXX This currently works also for emission classes and mixtures: specialize
        hmm = self.cmodel
        strout = "\nOverview of HMM:"
        strout += "\nNumber of states: " + str(hmm.N)
        strout += "\nNumber of mixture components: " + str(hmm.M)
        strout += "\nNumber of Output-Classes: " + str(hmm.cos)
        
        for k in range(hmm.N):
            print k
            state = ghmmwrapper.get_sstate(hmm, k)
            strout += "\n\nState number "+ str(k) + ":"
            strout += "\nInitial probability: " + str(state.pi)
            strout += "\n"+ str(hmm.M) + " density function(s):\n"
            weight = ""
            mue = ""
            u =  ""
            for outp in range(hmm.M):
                weight += str(ghmmwrapper.get_arrayd(state.c,outp)) + ", "
                mue += str(ghmmwrapper.get_arrayd(state.mue,outp)) + ", "
                u += str(ghmmwrapper.get_arrayd(state.u,outp)) + ", "
                strout += "  pdf component weights : " + str(weight) + "\n"
                strout += "  mean vector: " + str(mue) + "\n"
                strout += "  variance vector: " + str(u) + "\n"
                strout += "\nOutgoing transitions:"

                for i in range( state.out_states):
                    strout += "\ntransition to node " + str(ghmmwrapper.get_arrayint(state.out_id,i) )
                    for j in range(hmm.cos):
                        strout += "\n\tin class " + str(j) + " with probablity = " + str(ghmmwrapper.get_2d_arrayd(state.out_a,j,i))
                strout +=  "\nIngoing transitions:"
                for i in range(state.in_states):
                    strout += "\ntransition from node " + str(ghmmwrapper.get_arrayint(state.in_id,i) )
                    for j in range(hmm.cos):
                        strout += "\n\tin class "+str(j)+" with probablity = "+ str(ghmmwrapper.get_2d_arrayd(state.in_a,j,i))
            strout += "\nint fix:" + str(state.fix) + "\n"
        return strout




if __name__ == '__main__':
	m = HMMFromMatrices(DNA,DiscreteDistribution(DNA),
						[[0.0,1.0,0],[0.5,0.0,0.5],[1.0,0.0,0.0]],
						[[1.0,0.0,0.0,0.0],[0.0,0.5,0.5,0.0], [0.0,0.0,1.0,0.0]],
						[1.0,0,0])


    print m


	m2 = HMMFromMatrices(Float,GaussianDistribution(Float),
						[[0.0,1.0,0],[0.5,0.0,0.5],[1.0,0.0,0.0]],
						[[0.0,1.0],[-1.0,0.5], [1.0,0.2]],
						[1.0,0,0])


    print m2


	seq_c = ghmmwrapper.seq_d_read('test10.sqd')
	l = m2.loglikelihood_sqd(seq_c)

	for i in xrange(seq_c.seq_number):
		print ghmmwrapper.get_arrayd(l,i), 


	#m = HMMOpen("test.smo", modelIndex = 3) # Pick 3-rd model out of the smo fiel
	#m = HMMOpen("test.smo")

	#seqs = SequenceSetOpen('test.sqd')
	#l = m.baumWelch(seqs, 100, 0.001)
	#print l
