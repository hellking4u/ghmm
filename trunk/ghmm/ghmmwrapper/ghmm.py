#!/usr/bin/env python2.3
################################################################################
#
#       This file is part of the General Hidden Markov Model Library,
#       GHMM version __VERSION__, see http://ghmm.org
#
#       file:    ghmm.py
#       authors: Benjamin Georgi, Wasinee Rungsarityotin, Alexander Schliep,
#                Janne Grunau
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

   For technical reasons there can be two representations of an
   emission symbol: an external and an internal. The external
   representation is the view of the application using ghmm.py. The
   internal one is what is used in both ghmm.py and the ghmm
   C-library. Representations can coincide, but this is not
   guaranteed. Discrete alphabets of size k are represented as
   [0,1,2,...,k-1] internally.  It is the domain objects job to
   provide a mapping between representations in both directions.

   NOTE: Do not make assumptions about the internal
   representations. It might change.


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


4) The HMM: The HMM consists of two major components: A Markov chain
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

   If you want to store application specific data for each state you
   have to do it yourself.

   Subclasses of HMM implement specific types of HMM. The type depends
   on the EmissionDomain, the Distribution used, the specific
   extensions to the 'standard' HMMs and so forth


5) HMMFactory: This provides a way of constucting HMMs. Classes derived
   from HMMFactory allow to read HMMs from files, construct them explicitly
   from, for a discrete alphabet, transition matrix, emission matrix and prior
   or serve as the basis for GUI-based model building.

   There are several ways of using the HMMFactory.

   Static construction:

   HMMOpen(fileName) # Calls an object of type HMMOpen instantiated in ghmm
   HMMOpen(fileName, type=HMM.FILE_XML)
   HMMFromMatrices(emission_domain, distribution, A, B, pi) # B is a list of distribution parameters

   Examples:

   hmm = HMMOpen('some-hmm.xml')


   hmm.bla ....



"""

import ghmmwrapper
import ghmmhelper
import modhmmer
import re
import StringIO
import copy
import math
import sys
import os
import logging
from string import join
from textwrap import fill

# Initialize logging to stderr
#logging.basicConfig(format="%(asctime)s %(filename)s:%(lineno)d %(levelname)-5s - %(message)s")
log = logging.getLogger("GHMM")

# creating StreamHandler to stderr
hdlr = logging.StreamHandler(sys.stderr)

# setting message format
#fmt = logging.Formatter("%(name)s %(asctime)s %(filename)s:%(lineno)d %(levelname)s %(thread)-5s - %(message)s")
fmt = logging.Formatter("%(name)s %(filename)s:%(lineno)d  - %(message)s")
hdlr.setFormatter(fmt)

# adding handler to logger object
log.addHandler(hdlr)

# Set the minimal severity of a message to be shown. The levels in
# increasing severity are: DEBUG, INFO, WARNING, ERROR, CRITICAL

log.setLevel(logging.ERROR)
log.info( " I'm the ghmm in "+ __file__)

c_log = [log.critical, log.error, log.warning, log.info, log.debug]
def logwrapper(level, message):
    c_log[level](message)

ghmmwrapper.set_pylogging(logwrapper)

# Initialize global random number generator by system time
ghmmwrapper.ghmm_rng_init()
ghmmwrapper.time_seed()


#-------------------------------------------------------------------------------
#- Exceptions ------------------------------------------------------------------

class GHMMError(Exception):
    """Base class for exceptions in this module."""
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class UnknownInputType(GHMMError):
    def __init__(self,message):
       self.message = message
    def __str__(self):
        return repr(self.message)


class NoValidCDataType(GHMMError):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)


class badCPointer(GHMMError):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)


class SequenceCannotBeBuild(GHMMError):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class IndexOutOfBounds(GHMMError):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class InvalidModelParameters(GHMMError):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class GHMMOutOfDomain(GHMMError):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class UnsupportedFeature(GHMMError):
    def __init__(self,message):
       self.message = message
    def __str__(self):
        return repr(self.message)
    
class WrongFileType(GHMMError):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

#-------------------------------------------------------------------------------
#- constants -------------------------------------------------------------------
kNotSpecified            = ghmmwrapper.kNotSpecified
kLeftRight               = ghmmwrapper.kLeftRight
kSilentStates            = ghmmwrapper.kSilentStates
kTiedEmissions           = ghmmwrapper.kTiedEmissions
kHigherOrderEmissions    = ghmmwrapper.kHigherOrderEmissions
kBackgroundDistributions = ghmmwrapper.kBackgroundDistributions
kLabeledStates           = ghmmwrapper.kLabeledStates
kTransitionClasses       = ghmmwrapper.kTransitionClasses
kDiscreteHMM             = ghmmwrapper.kDiscreteHMM
kContinuousHMM           = ghmmwrapper.kContinuousHMM
kPairHMM                 = ghmmwrapper.kPairHMM
types = {
    kLeftRight:'kLeftRight',
    kSilentStates:'kSilentStates',
    kTiedEmissions:'kTiedEmissions',
    kHigherOrderEmissions:'kHigherOrderEmissions',
    kBackgroundDistributions:'kBackgroundDistributions',
    kLabeledStates:'kLabeledStates',
    kTransitionClasses:'kTransitionClasses',
    kDiscreteHMM:'kDiscreteHMM',
    kContinuousHMM:'kContinuousHMM',
    kPairHMM:'kPairHMM',
    }
#-------------------------------------------------------------------------------
#- EmissionDomain and derived  -------------------------------------------------
class EmissionDomain(object):
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
    def __init__(self, listOfCharacters, calphabet = None):
        """
        Creates an alphabet out of a listOfCharacters
        @param listOfCharacters: a list of strings (single characters most of
        the time), ints, or other objects that can be used as dictionary keys
        for a mapping of the external sequences to the internal representation
        Can alternatively be an C alphabet_s struct
        
        Note: Alphabets should be considered as imutable. That means the
        listOfCharacters and the mapping should never be touched after
        construction.
        """
        self.index = {} # Which index belongs to which character

        if calphabet is None:
            self.listOfCharacters = listOfCharacters
        else:
            self.listOfCharacters = [calphabet.getSymbol(i) for i in range(calphabet.size)]

        for i,c in enumerate(self.listOfCharacters):
            self.index[c] = i

        lens = {}
        try:
            for c in self.listOfCharacters:
                lens[len(c)] = 1
        except TypeError:
            self._lengthOfCharacters = None
        else:
            if len(lens) == 1:
                self._lengthOfCharacters = lens.keys()[0]
            else:
                self._lengthOfCharacters = None

        self.CDataType = "int" # flag indicating which C data type should be used


    def __str__(self):
        strout = ["<Alphabet:"]
        strout.append( str(self.listOfCharacters) +'>')

        return join(strout,'')
    
    def verboseStr(self):
        strout = ["GHMM Alphabet:\n"]
        strout.append("Number of symbols: " + str(len(self)) + "\n")
        strout.append("External: " + str(self.listOfCharacters) + "\n")
        strout.append("Internal: " + str(range(len(self))) + "\n")
        return join(strout,'')


    def __eq__(self,alph):
        if not isinstance(alph,Alphabet):
            return False
        else:    
            if self.listOfCharacters == alph.listOfCharacters and self.index == alph.index and self.CDataType==alph.CDataType:
                return True
            else:
                return False    

    def __len__(self):
        return len(self.listOfCharacters)

    def __hash__(self):
        #XXX rewrite 
        # defining hash and eq is not recommended for mutable types.
        # => listOfCharacters should be considered imutable
        return id(self)

    def size(self):
        #XXX
        """ Deprecated """
        log.warning( "Warning: The use of .size() is deprecated. Use len() instead.")
        return len(self.listOfCharacters)

        
    def internal(self, emission):
        """ Given a emission return the internal representation
        """
        return self.index[emission]


    def internalSequence(self, emissionSequence):
        """ Given a emission_sequence return the internal representation

            Raises KeyError
        """
        result = copy.deepcopy(emissionSequence)
        try:
            result = map(lambda i: self.index[i], result)
        except IndexError:
            raise KeyError
        return result


    def external(self, internal):
        """ Given an internal representation return the
            external representation

            Note: the internal code -1 always represents a gap character '-'

            Raises KeyError
        """
        if internal == -1:
            return "-"
        if internal < -1 or len(self.listOfCharacters) < internal:
            raise KeyError, "Internal symbol "+str(internal)+" not recognized."
        return self.listOfCharacters[internal]

    def externalSequence(self, internalSequence):
        """ Given a sequence with the internal representation return the
            external representation

            Raises KeyError
        """
        result = copy.deepcopy(internalSequence)
        try:
            result = map(lambda i: self.listOfCharacters[i], result)
        except IndexError:
            raise KeyError
        return result

    def isAdmissable(self, emission):
        """ Check whether emission is admissable (contained in) the domain
        """
        return emission in self.listOfCharacters

    def getExternalCharacterLength(self):
        """
        If all external characters are of the same length the length is
        returned. Otherwise None.
        @return: length of the external characters or None
        """
        return self._lengthOfCharacters

    def toCstruct(self):
        calphabet = ghmmwrapper.ghmm_alphabet(len(self), "<unused>")
        for i,symbol in enumerate(self.listOfCharacters):
            calphabet.setSymbol(i, str(symbol))

        return calphabet


DNA = Alphabet(['a','c','g','t'])
AminoAcids = Alphabet(['A','C','D','E','F','G','H','I','K','L',
                       'M','N','P','Q','R','S','T','V','W','Y'])
def IntegerRange(a,b):
    return Alphabet(range(a,b))


# To be used for labelled HMMs. We could use an Alphabet directly but this way it is more explicit.
class LabelDomain(Alphabet):    
    def __init__(self, listOfLabels, calphabet = None):
        Alphabet.__init__(self, listOfLabels, calphabet)


class Float(EmissionDomain):

    def __init__(self):
        self.CDataType = "double" # flag indicating which C data type should be used

    def __eq__(self, other):
        return isinstance(other, Float)

    def __hash__(self):
        # defining hash and eq is not recommended for mutable types.
        # for float it is fine because it is kind of state less
        return id(self)

    def isAdmissable(self, emission):
        """ Check whether emission is admissable (contained in) the domain
            raises GHMMOutOfDomain else
        """
        return isinstance(emission,float)



#-------------------------------------------------------------------------------
#- Distribution and derived  ---------------------------------------------------
class Distribution(object):
    """ Abstract base class for distribution over EmissionDomains
    """

    # add density, mass, cumuliative dist, quantils, sample, fit pars,
    # moments


class DiscreteDistribution(Distribution):
    """ A DiscreteDistribution over an Alphabet: The discrete distribution
        is parameterized by the vectors of probabilities.

    """
    def __init__(self, alphabet):
        self.alphabet = alphabet
        self.prob_vector = None

    def set(self, prob_vector):
        self.prob_vector = prob_vector

    def get(self):
        return self.prob_vector


class ContinuousDistribution(Distribution):
    pass

class UniformDistribution(ContinuousDistribution):
    def __init__(self, domain):
        self.emissionDomain = domain
        self.max = None
        self.min = None

    def set(self, (max, min)):
        self.max = max
        self.min = min

    def get(self):
        return (self.max, self.min)

class GaussianDistribution(ContinuousDistribution):
    # XXX attributes unused at this point
    def __init__(self, domain):
        self.emissionDomain = domain
        self.mu = None
        self.sigma = None

    def set(self, (mu, sigma)):
        self.mu = mu
        self.sigma = sigma

    def get(self):
        return (self.mu, self.sigma)

class TruncGaussianDistribution(GaussianDistribution):
    # XXX attributes unused at this point
    def __init__(self, domain):
        self.GaussianDistribution(self,domain)
        self.trunc = None

    def set(self, (mu, sigma,trunc)):
        self.mu = mu
        self.sigma = sigma
        self.trunc = trunc
        
    def get(self):
        return (self.mu, self.sigma, self.trunc) 

class GaussianMixtureDistribution(ContinuousDistribution):
    # XXX attributes unused at this point
    def __init__(self, domain):
        self.emissionDomain = domain
        self.M = None   # number of mixture components
        self.mu = []
        self.sigma = []
        self.weight = []

    def set(self, index, (mu, sigma,w)):
        pass

    def get(self):
        pass

class ContinuousMixtureDistribution(ContinuousDistribution):
    def __init__(self, domain):
        self.emissionDomain = domain
        self.M = 0   # number of mixture components
        self.components = []
        self.weight = []
        self.fix = []

    def add(self,w,fix,distribution):
        assert isinstance(distribution,ContinuousDistribution)
        self.M = self.M + 1
        self.weight.append(w)
        self.component.append(distribution)
        if isinstance(distribution,UniformDistribution): 
	# uniform distributions are fixed by definition
	  self.fix.append(1)            
	else:
          self.fix.append(fix)

    def set(self, index, w, fix, distribution):
        assert M > index
        assert isinstance(distribution,ContinuousDistribution)
        self.weight[i] = w
        self.components[i] = distribution
	if isinstance(distribution,UniformDistribution): 
	# uniform distributions are fixed by definition
	  self.fix[i](1)            
	else:
          self.fix[i](fix)

    def get(self,i):
        assert M > i
        return (self.weigth[i],self.fix[i],self.component[i])

    def check(self):        
        assert self.M == len(self.components)
        assert sum(self.weight) == 1
        assert sum(self.weight > 1) == 0
	assert sum(self.weight < 0) == 0        



#-------------------------------------------------------------------------------
#Sequence, SequenceSet and derived  ------------------------------------------

class EmissionSequence(object):
    """ An EmissionSequence contains the *internal* representation of
        a sequence of emissions. It also contains a reference to the
        domain where the emission orginated from.
    """

    def __init__(self, emissionDomain, sequenceInput, labelDomain = None, labelInput = None, ParentSequenceSet=None):

        self.emissionDomain = emissionDomain

        if ParentSequenceSet is not None:
            # optional reference to a parent SequenceSet. Is needed for reference counting
            #XXX exception
            assert isinstance(ParentSequenceSet,SequenceSet), "Error: Invalid reference. Only SequenceSet is valid."    
            self.ParentSequenceSet = ParentSequenceSet
        else:
            self.ParentSequenceSet = None

        if self.emissionDomain.CDataType == "int":
            # necessary C functions for accessing the ghmm_dseq struct
            self.sequenceAllocationFunction = ghmmwrapper.ghmm_dseq
            self.allocSingleSeq = ghmmwrapper.int_array_alloc
            self.seq_read = ghmmwrapper.ghmm_dseq_read
            self.seq_ptr_array_getitem = ghmmwrapper.dseq_ptr_array_getitem
            self.sequence_carray = ghmmhelper.list2int_array
        elif self.emissionDomain.CDataType == "double":
            # necessary C functions for accessing the ghmm_cseq struct
            self.sequenceAllocationFunction = ghmmwrapper.ghmm_cseq
            self.allocSingleSeq = ghmmwrapper.double_array_alloc
            self.seq_read = ghmmwrapper.ghmm_cseq_read
            self.seq_ptr_array_getitem = ghmmwrapper.cseq_ptr_array_getitem
            self.sequence_carray = ghmmhelper.list2double_array
        else:
            raise NoValidCDataType, "C data type " + str(self.emissionDomain.CDataType) + " invalid."


        # check if ghmm is build with asci sequence file support
        if isinstance(sequenceInput, str) or isinstance(sequenceInput, unicode):
            if ghmmwrapper.ASCI_SEQ_FILE:
                if  not os.path.exists(sequenceInput):
                     raise IOError, 'File ' + str(sequenceInput) + ' not found.'
                else:
                    i = ghmmwrapper.int_array_alloc(1)
                    tmp = self.seq_read(sequenceInput, i)
                    seq_number = ghmmwrapper.int_array_getitem(i, 0)
                    if seq_number > 0:
                        self.cseq = self.seq_ptr_array_getitem(tmp, 0)
                        for n in range(1, seq_number):
                            seq = self.seq_ptr_array_getitem(tmp, n)
                            del seq
                    else:
                        #XXX appropiate ecxception
                        raise IOError, 'File ' + str(sequenceInput) + ' not valid.'

                    ghmmwrapper.free(tmp)
                    ghmmwrapper.free(i)

            else:
                raise UnsupportedFeature("asci sequence files are deprecated. Please convert your files"
                                       + " to the new xml-format or rebuild the GHMM with"
                                       + " the conditional \"GHMM_OBSOLETE\".")

        #create a ghmm_dseq with state_labels, if the appropiate parameters are set
        elif isinstance(sequenceInput, list):
            internalInput = self.emissionDomain.internalSequence(sequenceInput)
            seq = self.sequence_carray(internalInput)
            self.cseq = self.sequenceAllocationFunction(seq, len(sequenceInput))

            if labelInput is not None and labelDomain is not None:
                assert len(sequenceInput)==len(labelInput), "Length of the sequence and labels don't match."
                assert isinstance(labelInput, list), "expected a list of labels."
                assert isinstance(labelDomain, LabelDomain), "labelDomain is not a LabelDomain class."
                
                self.labelDomain = labelDomain

                #translate the external labels in internal 
                internalLabel = self.labelDomain.internalSequence(labelInput)
                label = ghmmhelper.list2int_array(internalLabel)
                self.cseq.init_labels(label, len(internalInput))
                
        # internal use
        elif isinstance(sequenceInput, ghmmwrapper.ghmm_dseq) or isinstance(sequenceInput, ghmmwrapper.ghmm_cseq):
            if sequenceInput.seq_number > 1:
                raise badCPointer, "Use SequenceSet for multiple sequences."
            self.cseq = sequenceInput
            if labelDomain != None:
                self.labelDomain = labelDomain

        else:
            raise UnknownInputType, "inputType " + str(type(sequenceInput)) + " not recognized."


    def __del__(self):
        "Deallocation of C sequence struct."
        log.debug( "__del__ EmissionSequence " + str(self.cseq))

        # if a parent SequenceSet exits, we use cseq.subseq_free() to free memory
        if self.ParentSequenceSet is not None:
            self.cseq.subseq_free()


    def __len__(self):
        "Returns the length of the sequence."
        return self.cseq.getLength(0)


    def __setitem__(self, index, value):
        internalValue = self.emissionDomain.internal(value)
        self.cseq.setSymbol(0, index, internalValue)
		

    def __getitem__(self, index):
        """ Return the symbol at position 'index'. """
        if index < len(self):
            return self.cseq.getSymbol(0, index)
        else:
            raise IndexError    
    
    def getSeqLabel(self):
        if not ghmmwrapper.SEQ_LABEL_FIELD:
            raise UnsupportedFeature ("the seq_label field is obsolete. If you need it rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        return ghmmwrapper.long_array_getitem(self.cseq.seq_label,0)

    def setSeqLabel(self,value):
        if not ghmmwrapper.SEQ_LABEL_FIELD:
            raise UnsupportedFeature ("the seq_label field is obsolete. If you need it rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        ghmmwrapper.long_array_setitem(self.cseq.seq_label,0,value)
    
    def getStateLabel(self):
        """ Returns the labeling of the sequence in external representation"""
        if self.cseq.state_labels != None:
            iLabel = ghmmhelper.int_array2list(self.cseq.getLabels(0), self.cseq.getLabelsLength(0))
            return self.labelDomain.externalSequence(iLabel)
        else:
            #XXX appropiate exception
            raise IndexOutOfBounds(str(0) + " is out of bounds, only " + str(self.cseq.seq_number) + "labels")

    def getGeneratingStates(self):
        """
        Returns the state path from which the sequence was generated as a Python list.
        """
        l_state = []
        for j in range(ghmmwrapper.int_array_getitem(self.cseq.states_len,0) ):
            l_state.append(ghmmwrapper.int_matrix_getitem(self.cseq.states,0,j))

        return l_state

    def __str__(self):
        "Defines string representation."
        seq = self.cseq
        strout = []
       
        l = seq.getLength(0)
        if l <= 80:
       
            for j in range(l):
                strout.append(str( self.emissionDomain.external(self[j]) )   )
                if self.emissionDomain.CDataType == "double":
                    strout.append(" ")
        else:
            
            for j in range(0,5):
                strout.append(str( self.emissionDomain.external(self[j]) )   )
                if self.emissionDomain.CDataType == "double":
                    strout.append(" ")
            strout.append('...')
            for j in range(l-5,l):
                strout.append(str( self.emissionDomain.external(self[j]) )   )
                if self.emissionDomain.CDataType == "double":
                    strout.append(" ")

    	return join(strout,'')

    def verboseStr(self):
        "Defines string representation."
        seq = self.cseq
        strout = []
        strout.append("\nEmissionSequence Instance:\nlength " + str(seq.getLength(0)))
        strout.append(", weight " + str(seq.getWeight(0))  + ":\n")
        for j in range(seq.getLength(0)):
            strout.append(str( self.emissionDomain.external(self[j]) )   )
            if self.emissionDomain.CDataType == "double":
                strout.append(" ")

        # checking for labels
        if self.emissionDomain.CDataType == "int" and self.cseq.state_labels != None:            
            strout.append("\nState labels:\n")
            for j in range(seq.getLabelsLength(0)):
                strout.append(str( self.labelDomain.external(ghmmwrapper.int_matrix_getitem(seq.state_labels,0,j)))+ ", ")

    	return join(strout,'')


    def sequenceSet(self):
        """ Return a one-element SequenceSet with this sequence."""
        
        # in order to copy the sequence in 'self', we first create an empty SequenceSet and then
        # add 'self'
        seqSet = SequenceSet(self.emissionDomain, [])
        seqSet.cseq.add(self.cseq)
        return seqSet

    def write(self,fileName):
        "Writes the EmissionSequence into file 'fileName'."
        self.cseq.write(fileName)

    def setWeight(self, value):
        self.cseq.setWeight(0, value)
        self.cseq.total_w  = value
        
    def getWeight(self):
        return self.cseq.getWeight(0)

    def asSequenceSet(self):
        """returns a one element SequenceSet"""
        log.debug("EmissionSequence.asSequenceSet() -- begin " + repr(self.cseq))
        seq = self.sequenceAllocationFunction(1)

        # checking for state labels in the source C sequence struct
        if self.emissionDomain.CDataType == "int" and self.cseq.state_labels is not None:
            log.debug("EmissionSequence.asSequenceSet() -- found labels !")
            seq.calloc_state_labels()
            self.cseq.copyStateLabel(0, seq, 0)

        seq.setLength(0, self.cseq.getLength(0))
        seq.setSequence(0, self.cseq.getSequence(0))
        seq.setWeight(0, self.cseq.getWeight(0))
        
        # Above doesnt copy seq_id or seq_label or seq_w
        # XXX seq_id should be (long) int?
        seq_id = ghmmwrapper.double_array_getitem(self.cseq.seq_id, 0)
        ghmmwrapper.double_array_setitem(seq.seq_id, 0, seq_id)
        #if ghmmwrapper.SEQ_LABEL_FIELD:
        #    seq_label = ghmmwrapper.long_array_getitem(self.cseq.seq_label, i)
        #    ghmmwrapper.long_array_setitem(seq.seq_label, i, int(seq_label))

        log.debug("EmissionSequence.asSequenceSet() -- end " + repr(seq))
        return SequenceSetSubset(self.emissionDomain, seq, self)
        

class SequenceSet(object):
    def __init__(self, emissionDomain, sequenceSetInput, labelDomain = None, labelInput = None):
        self.emissionDomain = emissionDomain
        self.cseq = None

        if self.emissionDomain.CDataType == "int":
            # necessary C functions for accessing the ghmm_dseq struct
            self.sequenceAllocationFunction = ghmmwrapper.ghmm_dseq
            self.allocSingleSeq = ghmmwrapper.int_array_alloc
            self.seq_read = ghmmwrapper.ghmm_dseq_read
            self.seq_ptr_array_getitem = ghmmwrapper.dseq_ptr_array_getitem
            self.sequence_cmatrix = ghmmhelper.list2int_matrix
        elif self.emissionDomain.CDataType == "double":
            # necessary C functions for accessing the ghmm_cseq struct
            self.sequenceAllocationFunction = ghmmwrapper.ghmm_cseq
            self.allocSingleSeq = ghmmwrapper.double_array_alloc
            self.seq_read = ghmmwrapper.ghmm_cseq_read
            self.seq_ptr_array_getitem = ghmmwrapper.cseq_ptr_array_getitem
            self.sequence_cmatrix = ghmmhelper.list2double_matrix
        else:
            raise NoValidCDataType, "C data type " + str(self.emissionDomain.CDataType) + " invalid."


        # reads in the first sequence struct in the input file
        if isinstance(sequenceSetInput, str) or isinstance(sequenceSetInput, unicode):
            # check if ghmm is build with asci sequence file support
            if not ghmmwrapper.ASCI_SEQ_FILE:
                raise UnsupportedFeature ("asci sequence files are deprecated. Please convert your files to the new xml-format or rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
            else:
                if not os.path.exists(sequenceSetInput):
                    raise IOError, 'File ' + str(sequenceSetInput) + ' not found.'
                else:
                    i = ghmmwrapper.int_array_alloc(1)
                    tmp = self.seq_read(sequenceSetInput, i)
                    seq_number = ghmmwrapper.int_array_getitem(i, 0)
                    if seq_number > 0:
                        self.cseq = self.seq_ptr_array_getitem(tmp, 0)
                        for n in range(1, seq_number):
                            seq = self.seq_ptr_array_getitem(tmp, n)
                            del seq
                    else:
                        #XXX appropiate ecxception
                        raise IOError, 'File ' + str(sequenceSetInput) + ' not valid.'

                    ghmmwrapper.free(tmp)
                    ghmmwrapper.free(i)

        elif isinstance(sequenceSetInput, list):
            internalInput = [self.emissionDomain.internalSequence(seq) for seq in sequenceSetInput]
            (seq, lengths) = self.sequence_cmatrix(internalInput)
            lens = ghmmhelper.list2int_array(lengths)
                    
            self.cseq = self.sequenceAllocationFunction(seq, lens, len(sequenceSetInput))

            if isinstance(labelInput, list) and isinstance(labelDomain, LabelDomain): 
                assert len(sequenceSetInput)==len(labelInput), "no. of sequences and labels do not match."
                    
                self.labelDomain = labelDomain
                internalLabels = [self.labelDomain.internalSequence(oneLabel) for oneLabel in labelInput]
                (label,labellen) = ghmmhelper.list2int_matrix(internalLabels)
                lens = ghmmhelper.list2int_array(labellen)
                self.cseq.init_labels(label, lens)

        #internal use
        elif isinstance(sequenceSetInput, ghmmwrapper.ghmm_dseq) or isinstance(sequenceSetInput, ghmmwrapper.ghmm_cseq):
            log.debug("SequenceSet.__init__()", str(sequenceSetInput))
            self.cseq = sequenceSetInput
            if labelDomain is not None:
                self.labelDomain = labelDomain
                
        else:    
            raise UnknownInputType, "inputType " + str(type(sequenceSetInput)) + " not recognized."


    def __del__(self):
        "Deallocation of C sequence struct."
        log.debug( "__del__ SequenceSet " + str(self.cseq))


    def __str__(self):
        "Defines string representation."
        seq = self.cseq
        strout =  ["SequenceSet (N=" + str(seq.seq_number)+")"]

        
        if seq.seq_number <= 6:
            iter_list = range(seq.seq_number)
        else:
            iter_list = [0,2,'X',seq.seq_number-2,seq.seq_number-1]

        
        for i in iter_list:
            if i == 'X':
                strout.append('\n...\n')
            else:
                strout.append("\n  seq " + str(i)+ "(len=" + str(seq.getLength(i)) + ")\n")
                strout.append('    '+str(self[i]))


        return join(strout,'')


    def verboseStr(self):
        "Defines string representation."
        seq = self.cseq
        strout =  ["\nNumber of sequences: " + str(seq.seq_number)]

        for i in range(seq.seq_number):
            strout.append("\nSeq " + str(i)+ ", length " + str(seq.getLength(i)))
            strout.append(", weight " + str(seq.getWeight(i))  + ":\n")
            for j in range(seq.getLength(i)):
                if self.emissionDomain.CDataType == "int":
                    strout.append(str( self.emissionDomain.external(( ghmmwrapper.int_matrix_getitem(self.cseq.seq, i, j) )) ))
                elif self.emissionDomain.CDataType == "double":
                    strout.append(str( self.emissionDomain.external(( ghmmwrapper.double_matrix_getitem(self.cseq.seq, i, j) )) ) + " ")

            # checking for labels 
            if self.emissionDomain.CDataType == "int" and self.cseq.state_labels != None:            
                strout.append("\nState labels:\n")
                for j in range(seq.getLabelsLength(i)):
                    strout.append(str( self.labelDomain.external(ghmmwrapper.int_matrix_getitem(seq.state_labels,i,j))) +", ")

        return join(strout,'')


    def __len__(self):
        """ Return the number of sequences in the SequenceSet. """
        return self.cseq.seq_number

    def sequenceLength(self, i):
        """ Return the lenght of sequence 'i' in the SequenceSet """
        return self.cseq.getLength(i)

    def getWeight(self, i):
        """ Return the weight of sequence i. Weights are used in Baum-Welch"""
        return self.cseq.getWeight(i)

    def setWeight(self, i, w):
        """ Set the weight of sequence i. Weights are used in Baum-Welch"""
        ghmmwrapper.double_array_setitem(self.cseq.seq_w, i, w)
        
    def __getitem__(self, index):
        """ Return an EmissionSequence object initialized with a reference to 
        sequence 'index'.
        
        """
        # check the index for correct range
        if index >= self.cseq.seq_number:
            raise IndexError
        
        seq = self.cseq.get_singlesequence(index) 
        return EmissionSequence(self.emissionDomain, seq, ParentSequenceSet=self) 
    
    
    def getSeqLabel(self,index):
        if not ghmmwrapper.SEQ_LABEL_FIELD:
            raise UnsupportedFeature ("the seq_label field is obsolete. If you need it rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        return ghmmwrapper.long_array_getitem(self.cseq.seq_label,index)

    def setSeqLabel(self,index,value):
        if not ghmmwrapper.SEQ_LABEL_FIELD:
            raise UnsupportedFeature ("the seq_label field is obsolete. If you need it rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        ghmmwrapper.long_array_setitem(self.cseq.seq_label,index,value)

    def getGeneratingStates(self):
        """
        Returns the state paths from which the sequences were generated as a Python list of lists.
        """
        l_state = []
        for i in range(len(self)):
            ls_i = []
            for j in range(ghmmwrapper.int_array_getitem(self.cseq.states_len,i) ):
                ls_i.append(ghmmwrapper.int_matrix_getitem(self.cseq.states,i,j))
            l_state.append(ls_i)   

        return l_state


    def getSequence(self, index):
        """ Returns the index-th sequence in internal representation"""
        seq = []
        if self.cseq.seq_number > index:
            for j in range(self.cseq.getLength(index)):
                seq.append(self.cseq.getSymbol(index, j))
            return seq
        else:
            raise IndexOutOfBounds(str(index) + " is out of bounds, only " + str(self.cseq.seq_number) + "sequences")

    def getStateLabel(self,index):
        """ Returns the labeling of the index-th sequence in internal representation"""
        label = []
        if self.cseq.seq_number > index and self.cseq.state_labels != None:
            for j in range(self.cseq.getLabelsLength(index)):
                    label.append(self.labelDomain.external(ghmmwrapper.int_matrix_getitem(self.cseq.state_labels, index, j)))
            return label
        else:
            raise IndexOutOfBounds(str(0) + " is out of bounds, only " + str(self.cseq.seq_number) + "labels")
        

    def merge(self, emissionSequences): # Only allow EmissionSequence?
        """ 
             Merge 'emisisonSequences' with 'self'.
             'emisisonSequences' can either be an EmissionSequence or SequenceSet object.
        """

        if not isinstance(emissionSequences,EmissionSequence) and not isinstance(emissionSequences,SequenceSet):
            raise TypeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)
        
        self.cseq.add(emissionSequences.cseq)
        del(emissionSequences) # removing merged sequences

    def getSubset(self, seqIndixes):
        """ Returns a SequenceSet containing (references to) the sequences with the indixes in
            'seqIndixes'.
       
        """
        seqNumber = len(seqIndixes)
        seq = self.sequenceAllocationFunction(seqNumber)
        
        # checking for state labels in the source C sequence struct
        if self.emissionDomain.CDataType == "int" and self.cseq.state_labels is not None:
            
            log.debug( "SequenceSet: found labels !")
            seq.calloc_state_labels()

        for i in range(seqNumber):

            len_i = self.cseq.getLength(seqIndixes[i])
            
            seq.setSequence(i, self.cseq.getSequence(seqIndixes[i]))

            seq.setLength(i, len_i)

            # Above doesnt copy seq_id or seq_label or seq_w
            # XXX seq_id should be (long) int?
            seq_id = int(ghmmwrapper.double_array_getitem(self.cseq.seq_id, seqIndixes[i]))
            ghmmwrapper.double_array_setitem(seq.seq_id, i, seq_id)
            #if ghmmwrapper.SEQ_LABEL_FIELD:
            #    seq_label = ghmmwrapper.long_array_getitem(self.cseq.seq_label, i)
            #    ghmmwrapper.long_array_setitem(seq.seq_label, i, int(seq_label))

            seq.setWeight(i, self.cseq.getWeight(i))
             
            # setting labels if appropriate
            if self.emissionDomain.CDataType == "int" and self.cseq.state_labels is not None:
                self.cseq.copyStateLabel(seqIndixes[i], seq, seqIndixes[i])

        seq.seq_number = seqNumber
        
        return SequenceSetSubset(self.emissionDomain, seq, self)
        
    def write(self,fileName):
        "Writes (appends) the SequenceSet into file 'fileName'."
        self.cseq.write(fileName)

    def asSequenceSet(self):
        """conveinence function, returns only self"""
        return self

class SequenceSetSubset(SequenceSet):
    """ 
    SequenceSetSubset contains a subset of the sequences from a SequenceSet object.
    On the C side only the references are used.
    """
    def __init__(self, emissionDomain, sequenceSetInput, ParentSequenceSet , labelDomain = None, labelInput = None):
        # reference on the parent SequenceSet object
        log.debug("SequenceSetSubset.__init__ -- begin -", str(ParentSequenceSet))
        self.ParentSequenceSet = ParentSequenceSet
        SequenceSet.__init__(self, emissionDomain, sequenceSetInput, labelDomain, labelInput)

    def __del__(self):
        """ Since we do not want to deallocate the sequence memory, the destructor has to be
            overloaded.
        """
        log.debug( "__del__ SequenceSubSet " + str(self.cseq))
        
        if self.cseq is not None:
            self.cseq.subseq_free()

        # remove reference on parent SequenceSet object
        self.ParentSequenceSet = None



def SequenceSetOpen(emissionDomain, fileName):
    # XXX Name doof
    """ Reads a sequence file with multiple sequence sets. 

    Returns a list of SequenceSet objects.
    
    """

    if not os.path.exists(fileName):
        raise IOError, 'File ' + str(fileName) + ' not found.'

    if emissionDomain.CDataType == "int":
        readFile = ghmmwrapper.ghmm_dseq_read
        seqPtr   = ghmmwrapper.dseq_ptr_array_getitem
    elif emissionDomain.CDataType == "double":
        readFile = ghmmwrapper.ghmm_cseq_read
        seqPtr   = ghmmwrapper.cseq_ptr_array_getitem
    else:
        raise TypeError, "Invalid c data type " + str(emissionDomain.CDataType)

    # XXX Check whether swig can return setNr too        
    dArr = ghmmwrapper.int_array_alloc(1)

    structArray = readFile(fileName, dArr)
    setNr = ghmmwrapper.int_array_getitem(dArr, 0)

    # Add Unittest
    sequenceSets = [SequenceSet(emissionDomain, seqPtr(structArray, i)) for i in range(setNr)]
##    sequenceSets = []
##    for i in range(setNr):
##        seq = seqPtr(structArray,i)
##        sequenceSets.append(SequenceSet(emissionDomain, seq) )

##        # setting labels to NULL (XXX only for integer?. DONE by c-constructor) 
##        sequenceSets[i].cseq.state_labels = None
##        sequenceSets[i].cseq.state_labels_len = None
       
    ghmmwrapper.free(dArr)
    return  sequenceSets

def writeToFasta(seqSet,fn):
    """
    Writes a SequenceSet into a fasta file.
    """
    assert isinstance(seqSet, SequenceSet)
    f = open(fn,'w')
    
    for i in range(len(seqSet)):
        rseq = []
        for j in range(seqSet.sequenceLength(i)):
           rseq.append(str( seqSet.emissionDomain.external(( ghmmwrapper.int_matrix_getitem(seqSet.cseq.seq, i, j) )) ))
        
        f.write('>seq'+str(i)+'\n')
        f.write(fill(join(rseq,'') ))
        f.write('\n')

    f.close()        
   
    
    


#-------------------------------------------------------------------------------
# HMMFactory and derived  -----------------------------------------------------
class HMMFactory(object):
    """ A HMMFactory is the base class of HMM factories.
        A HMMFactory has just a constructor and a () method
    """


GHMM_FILETYPE_SMO = 'smo'
GHMM_FILETYPE_XML = 'xml'
GHMM_FILETYPE_HMMER = 'hmm'

# XXX Determine file type from file extension. Default unclear?

class HMMOpenFactory(HMMFactory):

    def __init__(self, defaultFileType=None):
        if defaultFileType:
            self.defaultFileType = defaultFileType

    def __call__(self, fileName, modelIndex = None):
        
        if not isinstance(fileName,StringIO.StringIO):
            if not os.path.exists(fileName):
                raise IOError, 'File ' + str(fileName) + ' not found.'
            
        # XML file: both new and old format
    	if self.defaultFileType == GHMM_FILETYPE_XML:
            # try to validate against ghmm.dtd
            if ghmmwrapper.ghmm_xmlfile_validate(fileName):
                return self.openNewXML(fileName, modelIndex)
            else:
                return self.openOldXML(fileName)
        elif self.defaultFileType == GHMM_FILETYPE_SMO:
            return self.openSMO(fileName, modelIndex)
        elif self.defaultFileType == GHMM_FILETYPE_HMMER:
            return self.openHMMER(fileName)
        else:
            raise TypeError, "Invalid file type " + str(self.defaultFileType)


    def openNewXML(self, fileName, modelIndex):
        # XXX Document me!
        # check the type of hmm
        # start the model

        file = ghmmwrapper.ghmm_xmlfile_parse(fileName)
        if file == None:
            log.debug( "XML has file format problems!")
            raise WrongFileType("file is not in GHMM xml format")

        nrModels = file.noModels
        modelType = file.modelType

        if (modelType & ghmmwrapper.kContinuousHMM):
            emission_domain = Float()
            distribution = ContinuousMixtureDistribution
            hmmClass = ContinuousMixtureHMM
            getPtr = ghmmwrapper.cmodel_ptr_array_getitem
            models = file.model.c

        elif ((modelType & ghmmwrapper.kDiscreteHMM)
              and not (modelType & ghmmwrapper.kTransitionClasses)
              and not (modelType & ghmmwrapper.kPairHMM)):
            emission_domain = 'd'
            distribution = DiscreteDistribution
            getPtr = ghmmwrapper.dmodel_ptr_array_getitem
            models = file.model.d
            if (modelType & ghmmwrapper.kLabeledStates):
                hmmClass = StateLabelHMM
            elif (modelType & ghmmwrapper.kDiscreteHMM):
                hmmClass = DiscreteEmissionHMM
            else:
                raise UnsupportedFeature, "Non-supported model type"

        else:
            raise UnsupportedFeature, "Non-supported model type"


        # read all models to list at first
        result = []
        for i in range(nrModels):
            cmodel = getPtr(models,i)
            if emission_domain is 'd': # XXX Uses first alphabet for all models in file
                emission_domain = Alphabet([], cmodel.alphabet)
                #print emission_domain
            if modelType & ghmmwrapper.kLabeledStates:
                labelDomain = LabelDomain([], cmodel.label_alphabet)
                #print labelDomain
                m = hmmClass(emission_domain, distribution(emission_domain), labelDomain, cmodel)
            else:
                m = hmmClass(emission_domain, distribution(emission_domain), cmodel)

            result.append(m)

        # for a single 
        if modelIndex != None:
            if modelIndex < nrModels:
                result = result[modelIndex]
            else:
                raise IndexOutOfBounds("the file %s has only %s models"% fileName, str(nrModels))
        elif nrModels == 1:
            result = result[0]


        del models
        return result

    def openOldXML(self, fileName):
        from ghmm_gato import xmlutil
        hmm_dom = xmlutil.HMM(fileName)
        emission_domain = hmm_dom.AlphabetType()

        if emission_domain == int:
            [alphabets, A, B, pi, state_orders] = hmm_dom.buildMatrices()

            emission_domain = Alphabet(alphabets)
            distribution = DiscreteDistribution(emission_domain)
            # build adjacency list

            # check for background distributions
            (background_dist, orders, code2name) = hmm_dom.getBackgroundDist()
            # (background_dist, orders) = hmm_dom.getBackgroundDist()
            bg_list = []
            # if background distribution exists, set background distribution here
            if background_dist != {}:
                # transformation to list for input into BackgroundDistribution,
                # ensure the rigth order
                for i in range(len(code2name.keys())-1):
                    bg_list.append(background_dist[code2name[i]])

                bg = BackgroundDistribution(emission_domain, bg_list)

            # check for state labels
            (label_list, labels) = hmm_dom.getLabels()
            if labels == ['None']:
                labeldom   = None
                label_list = None
            else:
                labeldom = LabelDomain(labels)

            m = HMMFromMatrices(emission_domain, distribution, A, B, pi, None, labeldom, label_list)

            # old xml is discrete, set appropiate flag
            m.cmodel.setModelTypeFlag(ghmmwrapper.kDiscreteHMM)

            if background_dist != {}:
                 ids = [-1]*m.N
                 for s in hmm_dom.state.values():
                      ids[s.index-1] = s.background # s.index ranges from [1, m.N]  

                 m.setBackground(bg, ids)
                 log.debug( "model_type %x" % m.cmodel.model_type)
                 log.debug("background_id" + str( ghmmhelper.int_array2list(m.cmodel.background_id, m.N)))
            else:
                 m.cmodel.bp = None
                 m.cmodel.background_id = None

            # check for tied states
            tied = hmm_dom.getTiedStates()
            if len(tied) > 0:
                m.setFlags(kTiedEmissions)
                m.cmodel.tied_to = ghmmhelper.list2int_array(tied)

            durations = hmm_dom.getStateDurations()
            if len(durations) == m.N: 
                log.debug("durations: " + str(durations))
                m.extendDurations(durations)

            return m

    def openSMO(self, fileName, modelIndex):
        # MO & SMO Files, format is deprecated
        # check if ghmm is build with smo support
        if not ghmmwrapper.SMO_FILE_SUPPORT:
            raise UnsupportedFeature ("smo files are deprecated. Please convert your files"
                                      "to the new xml-format or rebuild the GHMM with the"
                                      "conditional \"GHMM_OBSOLETE\".")
            
        (hmmClass, emission_domain, distribution) = self.determineHMMClass(fileName)

        log.debug("determineHMMClass = "+ str(  (hmmClass, emission_domain, distribution)))
        
        nrModelPtr = ghmmwrapper.int_array_alloc(1)
    	    
        # XXX broken since silent states are not supported by .smo file format 
        if hmmClass == DiscreteEmissionHMM:
            models = ghmmwrapper.ghmm_dmodel_read(fileName, nrModelPtr)
            getPtr = ghmmwrapper.dmodel_ptr_array_getitem
            base_model_type = ghmmwrapper.KDiscreteHMM
        else:
            models = ghmmwrapper.ghmm_cmodel_read(fileName, nrModelPtr)
            getPtr = ghmmwrapper.cmodel_ptr_array_getitem
            base_model_type = ghmmwrapper.kContinuousHMM

        nrModels = ghmmwrapper.int_array_getitem(nrModelPtr, 0)
        if modelIndex == None:
            result = []
            for i in range(nrModels):
                cmodel = getPtr(models, i)
                cmodel.setModelTypeFlag(base_model_type)
                m = hmmClass(emission_domain, distribution(emission_domain), cmodel)
                result.append(m)
        else:
            if modelIndex < nrModels:
                cmodel = getPtr(models, modelIndex)
                cmodel.setModelTypeFlag(base_model_type)
                result = hmmClass(emission_domain, distribution(emission_domain), cmodel)
            else:
                result = None

        return result

    def openHMMER(self, fileName):
        # HMMER format models
        h = modhmmer.hmmer(fileName)
        
        if h.m == 4:  # DNA model
            emission_domain = DNA
        elif h.m == 20:   # Peptide model    
            emission_domain = AminoAcids
        else:   # some other model
            emission_domain = IntegerRange(0,h.m)
        distribution = DiscreteDistribution(emission_domain)
        
        # XXX Probably slow for large matrices (Rewrite for 0.9)
        [A,B,pi,modelName] = h.getGHMMmatrices()
        return  HMMFromMatrices(emission_domain, distribution, A, B, pi, hmmName=modelName)

        
    def determineHMMClass(self, fileName):
        #
        # smo files. Obsolete 
        #
        file = open(fileName,'r')
        
        hmmRe = re.compile("^HMM\s*=")
        shmmRe = re.compile("^SHMM\s*=")
        mvalueRe = re.compile("M\s*=\s*([0-9]+)")
        densityvalueRe = re.compile("density\s*=\s*([0-9]+)")
        cosvalueRe = re.compile("cos\s*=\s*([0-9]+)")
        emission_domain = None

        while 1:
            l = file.readline()
            if not l:
                break
            l = l.strip()
            if len(l) > 0 and l[0] != '#': # Not a comment line
                hmm = hmmRe.search(l)
                shmm = shmmRe.search(l)
                mvalue = mvalueRe.search(l)
                densityvalue = densityvalueRe.search(l)
                cosvalue = cosvalueRe.search(l)
                
                if hmm != None:
                    if emission_domain != None and emission_domain != 'int':
                        log.error( "HMMOpenFactory:determineHMMClass: both HMM and SHMM? " + str(emission_domain))
                    else:
                        emission_domain = 'int'
                    
                if shmm != None:
                    if emission_domain != None and emission_domain != 'double':
                        log.error( "HMMOpenFactory:determineHMMClass: both HMM and SHMM? " + str(emission_domain))
                    else:
                        emission_domain = 'double'

                if mvalue != None:
                    M = int(mvalue.group(1))

                if densityvalue != None:
                    density = int(densityvalue.group(1))

                if cosvalue != None:
                    cos = int(cosvalue.group(1))

        file.close()
        if emission_domain == 'int':
            # only integer alphabet
            emission_domain = IntegerRange(0,M)
            distribution = DiscreteDistribution
            hmm_class = DiscreteEmissionHMM
            return (hmm_class, emission_domain, distribution)
        
        elif emission_domain == 'double':
            # M        number of mixture components
            # density  component type
            # cos      number of state transition classes
            if M == 1 and density == 0:
                emission_domain = Float()
                distribution = GaussianDistribution
                hmm_class = GaussianEmissionHMM            
                return (hmm_class, emission_domain, distribution)

            elif  M > 1 and density == 0:
                emission_domain = Float()
                distribution = GaussianMixtureDistribution
                hmm_class = GaussianMixtureHMM
                return (hmm_class, emission_domain, distribution)

            else:
                raise TypeError, "Model type can not be determined."

        return (None, None, None)
            
HMMOpenHMMER = HMMOpenFactory(GHMM_FILETYPE_HMMER) # read single HMMER model from file
HMMOpenSMO   = HMMOpenFactory(GHMM_FILETYPE_SMO)
HMMOpenXML   = HMMOpenFactory(GHMM_FILETYPE_XML)
# XXX Use HMMOpen for all, determine filetype from extension, filetype override for __call__ 
HMMOpen      = HMMOpenFactory(GHMM_FILETYPE_XML)


def readMultipleHMMERModels(fileName):
    """
        Reads a file containing multiple HMMs in HMMER format, returns list of
        HMM objects.

    """
    # XXX Integrate into HMMOPen, check for single hmm files
    if not os.path.exists(fileName):
        raise IOError, 'File ' + str(fileName) + ' not found.'
    
    modelList = []
    string = ""
    f = open(fileName,"r")
    
    res = re.compile("^//")
    stat = re.compile("^ACC\s+(\w+)")
    for line in f.readlines():
        string = string + line
        m = stat.match(line)
        if m:
            name = m.group(1)
            log.info( "Reading model " + str(name) + ".")
            
        match = res.match(line)
        if match:
            fileLike = StringIO.StringIO(string)
            modelList.append(HMMOpenHMMER(fileLike))
            string = ""
            match = None

    return modelList
    

class HMMFromMatricesFactory(HMMFactory):
    """ XXX Document matrix formats """

    #XXX: this should use the editing context
    def __call__(self, emissionDomain, distribution, A, B, pi, hmmName = None, labelDomain= None, labelList = None, densities = None):
        if isinstance(emissionDomain,Alphabet):
            
            assert emissionDomain == distribution.alphabet, "emissionDomain and distribution must be compatible"
            
            # checking matrix dimensions and argument validation, only some obvious errors are checked
            if not len(A) == len(A[0]):
                raise InvalidModelParameters, "A is not quadratic."
            if not len(pi) == len(A):
                raise InvalidModelParameters,  "Length of pi does not match length of A."
            if not len(A) == len(B):
                raise InvalidModelParameters, " Different number of entries in A and B."

            if (labelDomain is None and labelList is not None) or (labelList is None and labelList is not None):
                raise InvalidModelParameters, "Specify either both labelDomain and labelInput or neither."
            
            if isinstance(distribution,DiscreteDistribution):
                
                # HMM has discrete emissions over finite alphabet: DiscreteEmissionHMM

                # XXX Aufhuebschen
                cmodel = ghmmwrapper.ghmm_dmodel()
                cmodel.model_type = ghmmwrapper.kDiscreteHMM
                cmodel.N = len(A)
                cmodel.M = len(emissionDomain)
                cmodel.prior = -1 # No prior by default
                
                # tie groups are deactivated by default
                cmodel.tied_to = None
                
                # assign model identifier (if specified)
                if hmmName != None:
                    cmodel.name = hmmName
                else:
                    cmodel.name = ''

                states = ghmmwrapper.dstate_array_alloc(cmodel.N)

                silent_flag = 0
                silent_states = []

                tmpOrder = []

                #initialize states
                for i in range(cmodel.N):
                    state = ghmmwrapper.dstate_array_getRef(states, i)
                    # compute state order
                    if cmodel.M > 1:
                        order = math.log(len(B[i]), cmodel.M)-1
                    else:
                        order = len(B[i]) - 1
                        
                    log.debug( "order in state " + str(i) + " = " + str(order) )
                    # check or valid number of emission parameters
                    order = int(order)
                    if  cmodel.M**(order+1) == len(B[i]):
                        tmpOrder.append(order)
                    else:
                        raise InvalidModelParameters, "The number of "+str(len(B[i]))+ " emission parameters for state "+str(i)+" is invalid. State order can not be determined."
                    
                    state.b = ghmmhelper.list2double_array(B[i])
                                        
                    state.pi = pi[i]
                    
                    if (sum(B[i]) == 0 ): 
                        silent_states.append(1)
                        silent_flag = 4
                    else:
                        silent_states.append(0)

                    #set out probabilities
                    state.out_states, state.out_id, state.out_a = ghmmhelper.extract_out(A[i])

                    #set "in" probabilities
                    A_col_i = map( lambda x: x[i], A)
                    # Numarray use A[,:i]
                    state.in_states, state.in_id, state.in_a = ghmmhelper.extract_out(A_col_i)
                    #fix probabilities in reestimation, else 0
                    state.fix = 0

                cmodel.s = states
                if silent_flag == 4:
                    cmodel.model_type |= silent_flag
                    cmodel.silent = ghmmhelper.list2int_array(silent_states)

                cmodel.maxorder = max(tmpOrder)
                if cmodel.maxorder > 0:
                    log.debug( "Set kHigherOrderEmissions.")
                    cmodel.model_type |= kHigherOrderEmissions
                    cmodel.order = ghmmhelper.list2int_array(tmpOrder)

                # initialize lookup table for powers of the alphabet size,
                # speeds up models with higher order states
                powLookUp = [1] * (cmodel.maxorder+2)
                for i in range(1,len(powLookUp)):
                    powLookUp[i] = powLookUp[i-1] * cmodel.M
                cmodel.pow_lookup = ghmmhelper.list2int_array(powLookUp)
                
                # check for state labels
                if labelDomain is not None and labelList is not None:
                    if not isinstance(labelDomain,LabelDomain):
                        raise TypeError, "LabelDomain object required."
                    
                    cmodel.model_type |= kLabeledStates
                    m = StateLabelHMM(emissionDomain, distribution, labelDomain, cmodel)
                    m.setLabels(labelList)
                    return m
                else:    
                    return DiscreteEmissionHMM(emissionDomain, distribution, cmodel)

            else:
                raise GHMMError(type(distribution), "Not a valid distribution for Alphabet") 
        else:

            if isinstance(distribution,GaussianDistribution):
                
                cmodel = ghmmwrapper.ghmm_cmodel()
                cmodel.model_type = kContinuousHMM
                cmodel.M = 1 # Number of mixture componenent for emission distribution
                cmodel.prior = -1 # Unused

                # determining number of transition classes
                cos = 0
                if type(A[0][0]) == list:
                    cos = len(A)
                    cmodel.N = len(A[0])
                    
                    # allocating class switching context
                    cmodel.class_change_alloc()
                else: 
                    cos = 1
                    cmodel.N = len(A)
                    A = [A]
                
                cmodel.cos = ghmmhelper.classNumber(A)  # number of transition classes in GHMM

                log.debug( "cmodel.cos = "+str(cmodel.cos))

                states = ghmmwrapper.cstate_array_alloc(cmodel.N)

                # Switching functions and transition classes are handled
                # elswhere

                #initialize states
                for i in range(cmodel.N):

                    state = ghmmwrapper.cstate_array_getRef(states, i)
                    state.pi = pi[i]
                    state.M = 1

                    # allocate arrays of emmission parameters
                    state.c = ghmmhelper.list2double_array([1.0]) # Mixture weights. Unused
                    (mu, sigma) = B[i]
                    state.mue = ghmmhelper.list2double_array([mu]) #mu = mue in GHMM C-lib.
                    state.u = ghmmhelper.list2double_array([sigma])
                    state.c = ghmmhelper.list2double_array([1.0])
                    state.a = ghmmhelper.list2double_array([0.0])
                    
                    # mixture fixing deactivated by default
                    state.mixture_fix = ghmmhelper.list2int_array([0])

                    # setting densities types (all normal by default)                    
                    densities = ghmmwrapper.density_array_alloc(1)
                    state.density = densities
                    state.setDensity(0,0)                                
                    
                    #set out probabilities
                    trans = ghmmhelper.extract_out_cos(A, cmodel.cos, i) 
                    state.out_states = trans[0]
                    state.out_id = trans[1]
                    state.out_a = trans[2] 

                    #set "in" probabilities
                    trans = ghmmhelper.extract_in_cos(A,cmodel.cos, i) 
                    state.in_states = trans[0]
                    state.in_id = trans[1]
                    state.in_a = trans[2]
                
                    state.fix = 0 # if fix = 1 exclude the state probabilities from reestimation

                #append states to model
                cmodel.s = states
                m = GaussianEmissionHMM(emissionDomain, distribution, cmodel)
                m.cmodel.class_change = None
                return m
            
            if isinstance(distribution, GaussianMixtureDistribution):
                log.debug( "*** mixture model")
                
                # Interpretation of B matrix for the mixture case (Example with three states and two components each):
                #  B = [ 
                #      [ ["mu11","mu12"],["sig11","sig12"],["w11","w12"]   ],
                #      [  ["mu21","mu22"],["sig21","sig22"],["w21","w22"]  ],
                #      [  ["mu31","mu32"],["sig31","sig32"],["w31","w32"]  ],
                #      ]
                
                cmodel = ghmmwrapper.ghmm_cmodel()
                cmodel.M = len(B[0][0]) # Number of mixture componenents for emission distribution
                cmodel.prior = -1 # Unused
                
                # determining number of transition classes
                cos = 0
                if type(A[0][0]) == list:
                    cos = len(A)
                    cmodel.N = len(A[0])
                    # allocating class switching context
                    cmodel.class_change_alloc()
                    
                else: 
                    cos = 1
                    cmodel.N = len(A)
                    A = [A]
                    
                cmodel.cos = ghmmhelper.classNumber(A)  # number of transition classes in GHMM
                
                states = ghmmwrapper.cstate_array_alloc(cmodel.N)

                # Switching functions and transition classes are handled
                # elswhere

                #initialize states
                for i in range(cmodel.N):
                    state = ghmmwrapper.cstate_array_getRef(states,i)
                    state.pi = pi[i]
		    state.M = len(B[0][0])

                    # allocate arrays of emmission parameters
                    state.c = ghmmhelper.list2double_array([1.0]) # Mixture weights. Unused
                    mu_list = B[i][0]
                    sigma_list = B[i][1]
                    weight_list = B[i][2]
                    
                    state.mue = ghmmhelper.list2double_array(mu_list) #mu = mue in GHMM C-lib.
                    state.u = ghmmhelper.list2double_array(sigma_list)
                    state.c = ghmmhelper.list2double_array(weight_list)
		    state.a = ghmmhelper.list2double_array([0.0] * state.M)

                    # setting densities types (all normal by default)                    
                    densities = ghmmwrapper.density_array_alloc(cmodel.M)
                    state.density = densities
                    for j in range(state.M):
                        state.setDensity(j,0)

                    # mixture fixing deactivated by default
                    state.mixture_fix = ghmmhelper.list2int_array([0] * state.M)
                    
                    #set out probabilities
                    trans = ghmmhelper.extract_out_cos(A, cmodel.cos, i) 
                    state.out_states = trans[0]
                    state.out_id = trans[1]
                    state.out_a = trans[2] 

                    #set "in" probabilities
                    trans = ghmmhelper.extract_in_cos(A,cmodel.cos, i) 
                    state.in_states = trans[0]
                    state.in_id = trans[1]
                    state.in_a = trans[2]
                
                    state.fix = 0 # if fix = 1, exclude state's probabilities from reestimation                

                #append states to model
                cmodel.s = states
		
		
                return GaussianMixtureHMM(emissionDomain, distribution, cmodel)
                
            elif isinstance(distribution, ContinuousMixtureDistribution):
                log.debug( "*** general mixture model")
                 
                # Interpretation of B matrix for the mixture case (Example with three states and two components each):
                #  B = [ 
                #      [ ["mu11","mu12"],["sig11","sig12"],["a11","a12"],["w11","w12"]   ],
                #      [  ["mu21","mu22"],["sig21","sig22"],["w21","w22"]  ],
                #      [  ["mu31","mu32"],["sig31","sig32"],["w31","w32"]  ],
                #      ]

                assert densities != None, "Continuous Mixture Distributions need a density type array"
                 
                cmodel = ghmmwrapper.ghmm_cmodel()
                cmodel.M = len(B[0][0]) # Number of mixture componenents for emission distribution
                cmodel.prior = 0 # Unused
                
                # determining number of transition classes
                cos = 0
                if type(A[0][0]) == list:
                    cos = len(A)
                    cmodel.N = len(A[0])
                    # allocating class switching context
                    cmodel.class_change_alloc()
                    
                else: 
                    cos = 1
                    cmodel.N = len(A)
                    A = [A]
                    
                cmodel.cos = ghmmhelper.classNumber(A)  # number of transition classes in GHMM
                
                states = ghmmwrapper.cstate_array_alloc(cmodel.N)

                # Switching functions and transition classes are handled
                # elswhere

                #initialize states
                for i in range(cmodel.N):
                    state = ghmmwrapper.cstate_array_getRef(states,i)
                    state.pi = pi[i]
		    state.M = len(B[0][0])

                    # allocate arrays of emmission parameters
                    state.c = ghmmhelper.list2double_array([1.0]) # Mixture weight
                    mu_list = B[i][0]
                    sigma_list = B[i][1]
                    a_list = B[i][2]
                    weight_list = B[i][3]
                    
                    state.mue = ghmmhelper.list2double_array(mu_list) #mu = mue in GHMM C-lib.
                    state.u = ghmmhelper.list2double_array(sigma_list)
                    state.c = ghmmhelper.list2double_array(weight_list)
		    state.a = ghmmhelper.list2double_array(a_list)

                    # setting densities types (all normal by default)                    
                    densit = ghmmwrapper.density_array_alloc(cmodel.M)
                    state.density = densit

                    mix_fix = [0] * state.M

                    for j in range(state.M):                   
                      state.setDensity(j,densities[i][j])
                      if densities[i][j] == ghmmwrapper.uniform:
                        mix_fix[j] = 1

                    # mixture fixing deactivated by default
                    state.mixture_fix = ghmmhelper.list2int_array(mix_fix)
                    
                    #set out probabilities
                    trans = ghmmhelper.extract_out_cos(A, cmodel.cos, i) 
                    state.out_states = trans[0]
                    state.out_id = trans[1]
                    state.out_a = trans[2] 

                    #set "in" probabilities
                    trans = ghmmhelper.extract_in_cos(A,cmodel.cos, i) 
                    state.in_states = trans[0]
                    state.in_id = trans[1]
                    state.in_a = trans[2]
                
                    state.fix = 0 # if fix = 1, exclude state's probabilities from reestimation                

                #append states to model
                cmodel.s = states
		
		
                return ContinuousMixtureHMM(emissionDomain, distribution, cmodel)
            else:
                raise GHMMError(type(distribution),
                                "Cannot construct model for this domain/distribution combination") 


HMMFromMatrices = HMMFromMatricesFactory()

#-------------------------------------------------------------------------------
#- Background distribution

class BackgroundDistribution(object):
    """ XXX doc string """
    def __init__(self, emissionDomain, bgInput):
        
        if type(bgInput) == list:
            self.emissionDomain = emissionDomain
            distNum = len(bgInput)
        
            order = ghmmwrapper.int_array_alloc(distNum)
            b = ghmmwrapper.double_matrix_alloc_row(distNum)
            for i in range(distNum):
                if len(emissionDomain) > 1:
                    o = math.log(len(bgInput[i]), len(emissionDomain)) - 1
                else:
                    o = len(bgInput[i]) - 1
                         
                assert (o % 1) == 0, "Invalid order of distribution " + str(i) + ": " + str(o)

                ghmmwrapper.int_array_setitem(order, i, int(o))
                # dynamic allocation, rows have different lenghts
                b_i = ghmmhelper.list2double_array(bgInput[i])
                ghmmwrapper.double_matrix_set_col(b, i, b_i)
    
            self.cbackground = ghmmwrapper.ghmm_dbackground(distNum, len(emissionDomain), order, b)

        elif isinstance(bgInput, ghmmwrapper.background_distributions):
            self.cbackground = bgInput
            self.emissionDomain = emissionDomain
             
        else:
            raise TypeError, "Input type "+str(type(bgInput)) +" not recognized."    

    def __del__(self):
        log.debug( "__del__ BackgroundDistribution " + str(self.cbackground))
        del self.cbackground
        self.cbackground = None
    
    def __str__(self):
        outstr = 'BackgroundDistribution (N= '+str(self.cbackground.n)+'):\n'
        outstr += str(self.emissionDomain) + "\n"
        d = ghmmhelper.double_matrix2list(self.cbackground.b, self.cbackground.n, len(self.emissionDomain))
        outstr += "Distributions:\n"   
        f = lambda x: "%.2f" % (x,)  # float rounding function 
        
        for i in range(self.cbackground.n):
            outstr += '  '+str(i+1) + ":(order= " + str(self.cbackground.getOrder(i))+"): "
            outstr += " "+join(map(f,d[i]),', ')+"\n"
        return outstr


    def verboseStr(self):
        outstr = "BackgroundDistribution instance:\n"
        outstr += "Number of distributions: " + str(self.cbackground.n)+"\n\n"
        outstr += str(self.emissionDomain) + "\n"
        d = ghmmhelper.double_matrix2list(self.cbackground.b, self.cbackground.n, len(self.emissionDomain))
        outstr += "Distributions:\n"   
        for i in range(self.cbackground.n):
            outstr += "  Order: " + str(self.cbackground.getOrder(i))+"\n"
            outstr += "  " + str(i+1) +": "+str(d[i])+"\n"
        return outstr
        
    def tolist(self):
        dim = self.cbackground.m
        distNum = self.cbackground.n
        orders = ghmmhelper.int_array2list(self.cbackground.order, distNum)
        B = []
        for i in xrange(distNum):
             order = orders[i]
             size = int(pow(m,(order+1)))
             b = [0.0]*size
             for j in xrange(size):
                  b[j] = ghmmwrapper.double_matrix_getitem(self.cbackground.b,i,j)
             B.append(b)
        return (distNum,orders,B)

#-------------------------------------------------------------------------------
#- HMM and derived  
class HMM(object):
    """ The HMM base class. 
        All functions where the C signatures allows it will 
        be  defined in here. Unfortunately there stil is a lot of overloading going on in 
        derived classes.
        
        Generic features (these apply to all derived classes):
            - Forward algorithm
            - Viterbi algorithm
            - Baum-Welch training
            - HMM distance metric
            - ... 
    
    """
    def __init__(self, emissionDomain, distribution, cmodel):
        self.emissionDomain = emissionDomain
        self.distribution = distribution
        self.cmodel = cmodel
        
        self.N = self.cmodel.N  # number of states
        self.M = self.cmodel.M  # number of symbols / mixture components


    def __del__(self):
        """ Deallocation routine for the underlying C data structures. """
        log.debug( "__del__ HMM" + str(self.cmodel))
        del self.cmodel
        self.cmodel = None

    def loglikelihood(self, emissionSequences):
        """ Compute log( P[emissionSequences| model]) using the forward algorithm
            assuming independence of the sequences in emissionSequences

            emissionSequences can either be a SequenceSet or a Sequence

            Result: log( P[emissionSequences| model]) of type float which is
            computed as \sum_{s} log( P[s| model]) when emissionSequences
            is a SequenceSet

        Note: The implementation does not compute the full forward matrix since we are only interested
              in the likelihoods in this case.
        """        
        return sum(self.loglikelihoods(emissionSequences))


    def loglikelihoods(self, emissionSequences): 
        """ Compute a vector ( log( P[s| model]) )_{s} of log-likelihoods of the
            individual emission_sequences using the forward algorithm

            emission_sequences is of type SequenceSet

            Result: log( P[emissionSequences| model]) of type float
                    (numarray) vector of floats

        """
        log.debug("HMM.loglikelihoods() -- begin")
        emissionSequences = emissionSequences.asSequenceSet()
        seqNumber = len(emissionSequences)

        likelihood = ghmmwrapper.double_array_alloc(1)
        likelihoodList = []

        for i in range(seqNumber):
            seq = emissionSequences.cseq.getSequence(i)
            tmp = emissionSequences.cseq.getLength(i)
            
            ret_val = self.cmodel.logp(seq, tmp, likelihood)
            if ret_val == -1:
                
                log.warning("forward returned -1: Sequence"+str(i)+"cannot be build.")
                # XXX Eventually this should trickle down to C-level
                # Returning -DBL_MIN instead of infinity is stupid, since the latter allows
                # to continue further computations with that inf, which causes
                # things to blow up later.
                # cmodel.logp() could do without a return value if -Inf is returned
                # What should be the semantics in case of computing the likelihood of
                # a set of sequences
                likelihoodList.append(-float('Inf'))
            else:
                likelihoodList.append(ghmmwrapper.double_array_getitem(likelihood,0))

        ghmmwrapper.free(likelihood)
        del emissionSequences
        log.debug("HMM.loglikelihoods() -- end")
        return likelihoodList

    ## Further Marginals ...

    def pathPosterior(self, sequence, path):
        """ Returns the log posterior probability for 'path' having generated 'sequence'.
            
            CAVEAT: statePosterior needs to calculate the complete forward and
            backward matrices. If you are interested in multiple paths it would
            be more efficient to use the 'posterior' function directly and not
            multiple calls to pathPosterior
        """
        # implemented in derived classes 
        raise NotImplementedError
        

    def statePosterior(self, sequence, state, time):
        """ Return the log posterior probability for being at 'state' at time 'time' in 'sequence'.
            
            CAVEAT: statePosterior needs to calculate the complete forward and backward matrices. If 
            you are interested in multiple states it would be more efficient to use the posterior function
            directly and not multiple calls to statePosterior
        
        """
        # implemented in derived classes 
        raise NotImplementedError


    def posterior(self, sequence):
        """ Posterior distribution matrix for 'sequence'.
        
        """
        # implemented in derived classes 
        raise NotImplementedError



    def logprob(self, emissionSequence, stateSequence):
        """log P[ emissionSequence, stateSequence| m] 
        """
        # implemented in derived classes.
        raise NotImplementedError

    # The functions for model training are defined in the derived classes.
    def baumWelch(self, trainingSequences, nrSteps, loglikelihoodCutoff):
        raise NotImplementedError

    def baumWelchSetup(self, trainingSequences, nrSteps):
        raise NotImplementedError

    def baumWelchStep(self, nrSteps, loglikelihoodCutoff):
        raise NotImplementedError
    
    def baumWelchDelete(self):
        raise NotImplementedError
        
    # extern double ghmm_c_prob_distance(smodel *cm0, smodel *cm, int maxT, int symmetric, int verbose);
    def distance(self, model, seqLength):
        """ Returns the distance between 'self.cmodel' and 'model'.   """
        return self.cmodel.prob_distance(model.cmodel, seqLength, 0, 0)


    def forward(self, emissionSequence):
        """
            Result: the (N x T)-matrix containing the forward-variables
                    and the scaling vector
        """
        log.debug("HMM.forward -- begin")
        # XXX Allocations should be in try, except, finally blocks
        # to assure deallocation even in the case of errrors.
        # This will leak otherwise.
        seq = emissionSequence.cseq.getSequence(0)        

        unused = ghmmwrapper.double_array_alloc(1) # Dummy return value for forward()    	
     
        t = len(emissionSequence)
        calpha = ghmmwrapper.double_matrix_alloc(t, self.N)
        cscale = ghmmwrapper.double_array_alloc(t)

        error = self.cmodel.forward(seq, t, calpha, cscale, unused)
        if error == -1:
            log.error( "forward finished with -1: EmissionSequence cannot be build.")

        # translate alpha / scale to python lists 
        pyscale = ghmmhelper.double_array2list(cscale, t)
        pyalpha = ghmmhelper.double_matrix2list(calpha, t, self.N)
        
        # deallocation
        ghmmwrapper.free(unused)
        ghmmwrapper.free(cscale)
        ghmmwrapper.double_matrix_free(calpha, t)

        log.debug("HMM.forward -- end")
        return pyalpha, pyscale
        

    def backward(self, emissionSequence, scalingVector):
        """
            Result: the (N x T)-matrix containing the backward-variables
        """
        log.debug("HMM.backward -- begin")
        seq = emissionSequence.cseq.getSequence(0)
        
        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2double_array(scalingVector)
        
        # alllocating beta matrix
        t = len(emissionSequence)
        cbeta = ghmmwrapper.double_matrix_alloc(t, self.N)
        
        error = self.cmodel.backward(seq,t,cbeta,cscale)
        if error == -1:
            log.error( "backward finished with -1: EmissionSequence cannot be build.")
            
        pybeta = ghmmhelper.double_matrix2list(cbeta,t,self.N)

        # deallocation
        ghmmwrapper.free(cscale)
        ghmmwrapper.double_matrix_free(cbeta,t)

        log.debug("HMM.backward -- end")
        return pybeta


    def viterbi(self, eseqs):
        """ Compute the Viterbi-path for each sequence in emissionSequences

            emission_sequences can either be a SequenceSet or an EmissionSequence

            Result: [q_0, ..., q_T] the viterbi-path of emission_sequences is an emmissionSequence
            object, [[q_0^0, ..., q_T^0], ..., [q_0^k, ..., q_T^k]} for a k-sequence
                    SequenceSet
        """
        log.debug("HMM.viterbi() -- begin")
        emissionSequences = eseqs.asSequenceSet()

        seqNumber = len(emissionSequences)        

        log_p = ghmmwrapper.double_array_alloc(1)

        allLogs = []
        allPaths = []
        for i in range(seqNumber):
            seq = emissionSequences.cseq.getSequence(i)
            seq_len = emissionSequences.cseq.getLength(i)
            
            if seq_len > 0:
                viterbiPath = self.cmodel.viterbi(seq,seq_len,log_p)
            else:
                viterbiPath = None

            onePath = []
            
            # for model types without possible silent states
            # the length of the viterbi path is known
            if not self.hasFlags(kSilentStates):
                onePath = ghmmhelper.int_array2list(viterbiPath, seq_len)
            
            # in the silent case we have to append as long as the path is
            # positive because the final path position is marked with a -1
            # in the following array entry.
            else:
                for j in range(seq_len * self.N): # maximum length of a viterbi path for a silent model
                    d = ghmmwrapper.int_array_getitem(viterbiPath, j)
                    if d >= 0:
                        onePath.append(d)
                    else:
                        break

            allPaths.append(onePath)
            allLogs.append(ghmmwrapper.double_array_getitem(log_p, 0))
            ghmmwrapper.free(viterbiPath) 

        ghmmwrapper.free(log_p)

        log.debug("HMM.viterbi() -- end")
        if seqNumber > 1:
            return allPaths, allLogs
        else:
            return allPaths[0], allLogs[0]


    def sample(self, seqNr ,T, seed = 0):
        """ Sample emission sequences 
                seqNr = number of sequences to be sampled
                T = length of each sequence
                seed = initialization value for rng, default 0 means 

        """
        seqPtr = self.cmodel.generate_sequences(seed,T,seqNr,self.N)
        return SequenceSet(self.emissionDomain,seqPtr)
        

    def sampleSingle(self, T, seed = 0):
        """ Sample a single emission sequence of length at most T.
            Returns a Sequence object.
        """
        log.debug("HMM.sampleSingle() -- begin")
        seqPtr = self.cmodel.generate_sequences(seed,T,1,self.N)
        log.debug("HMM.sampleSingle() -- end")
        return EmissionSequence(self.emissionDomain,seqPtr)

    def clearFlags(self, flags):
        """ Clears one or more model type flags. Use with care.
        """
        log.debug("clearFlags: " + self.printtypes(flags))
        self.cmodel.model_type &= ~flags

    def hasFlags(self, flags):
        """ Checks if the model has one or more model type flags set
        """
        return self.cmodel.model_type & flags

    def setFlags(self, flags):
        """ Sets one or more model type flags. Use with care.
        """
        log.debug("setFlags: " + self.printtypes(flags))
        self.cmodel.model_type |= flags

    def state(self, stateLabel):
        """ Given a stateLabel return the integer index to the state 

        """
        raise NotImplementedError

    def getInitial(self, i):
        """ Accessor function for the initial probability \pi_i """
        state = self.cmodel.getState(i)
        return state.pi

    def setInitial(self, i, prob, fixProb=False):
        """ Accessor function for the initial probability \pi_i
            If 'fixProb' = True \pi will be rescaled to 1 with 'pi[i]' fixed to the
            arguement value of 'prob'.

        """
        state = self.cmodel.getState(i)
        old_pi = state.pi
        state.pi = prob

        # renormalizing pi, pi(i) is fixed on value 'prob'
        if fixProb:
            coeff = (1.0 - old_pi) / prob
            for j in range(self.N):
                if i != j:
                    state = self.cmodel.getState(j)
                    p = state.pi
                    state.pi = p / coeff

    def getTransition(self, i, j):
        """ Accessor function for the transition a_ij """
        state = self.cmodel.getState(i)	

        # ensure proper indices
        # XXX IndexError
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."
        assert 0 <= j < self.N, "Index " + str(j) + " out of bounds."

        transition = 0.0
        for i in xrange(state.out_states):
            stateId = state.getOutState(i)
            if stateId == j:
                transition = state.getOutProb(i)
                break
        return transition

    def setTransition(self, i, j, prob):
        """ Accessor function for the transition a_ij. """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."
        assert 0 <= j < self.N, "Index " + str(j) + " out of bounds."

        # XXX Need to check that (i,j) is a transition, return IndexError else

        self.cmodel.set_transition(i, j, prob)


    def getEmission(self, i):
        """ Accessor function for the emission distribution parameters of state 'i'.

            For discrete models the distribution over the symbols is returned,
            for continuous models a matrix of the form
            [ [mu_1, sigma_1, weight_1] ... [mu_M, sigma_M, weight_M]  ] is returned.

        """
        raise NotImplementedError

    def setEmission(self, i, distributionParemters):
        """ Set the emission distribution parameters

            Defined in derived classes.
         """
        raise NotImplementedError

    # XXX asMatrices
    def toMatrices(self):
        "To be defined in derived classes."
        raise NotImplementedError


    def normalize(self):
        """ Normalize transition probs, emission probs (if applicable)

            Defined in derived classes.
        """
        raise NotImplementedError


    def randomize(self, noiseLevel):
        """ """
        raise NotImplementedError

    def write(self,fileName):
        """ Writes HMM to file 'fileName'.

        """
        self.cmodel.write_xml(fileName)


    def printtypes(self, model_type):
        strout = []
        if model_type == kNotSpecified:
            return 'kNotSpecified'
        for k in types.keys():
            if model_type & k:
                strout.append(types[k])
        return ' '.join(strout)


def HMMwriteList(fileName, hmmList, fileType=GHMM_FILETYPE_XML):
    if (fileType == GHMM_FILETYPE_XML):
        if os.path.exists(fileName):
            log.warning( "HMMwriteList: File " + str(fileName) + " already exists. Model will be overwritted.")
        models = ghmmwrapper.cmodel_ptr_array_alloc(len(hmmList))
        for i, model in enumerate(hmmList):
            ghmmwrapper.cmodel_ptr_array_setitem(models, i, model.cmodel)
        ghmmwrapper.ghmm_cmodel_xml_write(models, fileName, len(hmmList))
        ghmmwrapper.free(models)
    elif (fileType==GHMM_FILETYPE_SMO):
        raise WrongFileType("the smo file format is deprecated, use xml instead")
    else:
        raise WrongFileType("unknown file format" + str(fileType))
        
   
class DiscreteEmissionHMM(HMM):
    """ HMMs with discrete emissions.
        Optional features:
         - silent states
         - higher order states
         - parameter tying in training
         - background probabilities in training
    """

    def __init__(self, emissionDomain, distribution, cmodel):
        HMM.__init__(self, emissionDomain, distribution, cmodel)

        self.model_type = self.cmodel.model_type  # model type
        self.maxorder = self.cmodel.maxorder
        self.background = None


    def __str__(self):
        hmm = self.cmodel
        strout = [str(self.__class__.__name__)]
        if self.cmodel.name:
            strout.append( " " + str(self.cmodel.name))
        strout.append(  "(N= "+ str(hmm.N))
        strout.append(  ", M= "+ str(hmm.M)+')\n')
        
        f = lambda x: "%.2f" % (x,) # float rounding function 
        
        if self.hasFlags(kHigherOrderEmissions):
            order = ghmmhelper.int_array2list(self.cmodel.order, self.N)
        else:
            order = [0]*hmm.N
        
        if hmm.N <= 4:
            iter_list = range(self.N)
        else:
            iter_list = [0,1,'X',hmm.N-2,hmm.N-1]            

        for k in iter_list:
            if k == 'X':
                strout.append('\n  ...\n\n')
                continue
        
            state = hmm.getState(k)
            strout.append( "  state "+ str(k) +' (')
            if order[k] > 0:
                strout.append( 'order= '+ str(order[k])+',')


            strout.append( "initial= " + f(state.pi)+')\n')
            strout.append( "    Emissions: ")
            for outp in range(hmm.M**(order[k]+1)):
                strout.append(f(ghmmwrapper.double_array_getitem(state.b,outp)))
                if outp < hmm.M**(order[k]+1)-1:
                    strout.append( ', ')
                else:    
                    strout.append('\n')

            strout.append( "    Transitions:")
            #trans = [0.0] * hmm.N
            for i in range( state.out_states):
                strout.append( " ->" + str( state.getOutState(i))+' ('+ f(ghmmwrapper.double_array_getitem(state.out_a,i) ) +')' )
                if i < state.out_states-1:
                    strout.append( ',')
                #strout.append(" with probability " + str(ghmmwrapper.double_array_getitem(state.out_a,i)))

            strout.append('\n')

        return join(strout,'')


    
    def verboseStr(self):
        hmm = self.cmodel
        strout = ["\nGHMM Model\n"]
        strout.append( "Name: " + str(self.cmodel.name))
        strout.append( "\nModelflags: "+ self.printtypes(self.cmodel.model_type))
        strout.append(  "\nNumber of states: "+ str(hmm.N))
        strout.append(  "\nSize of Alphabet: "+ str(hmm.M))
        if self.hasFlags(kHigherOrderEmissions):
            order = ghmmhelper.int_array2list(self.cmodel.order, self.N)
        else:
            order = [0]*hmm.N

        for k in range(hmm.N):
            state = hmm.getState(k)
            strout.append( "\n\nState number "+ str(k) +":")
            strout.append( "\nState order: " + str(order[k]))
            strout.append( "\nInitial probability: " + str(state.pi))
            #strout.append("\nsilent state: " + str(self.cmodel.silent[k]))
            strout.append( "\nOutput probabilites: ")
            for outp in range(hmm.M**(order[k]+1)):
                strout.append(str(ghmmwrapper.double_array_getitem(state.b,outp)))
                if outp % hmm.M == hmm.M-1:
                    strout.append( "\n")
                else:
                    strout.append( ", ")

            strout.append( "\nOutgoing transitions:")
            for i in range( state.out_states):
                strout.append( "\ntransition to state " + str( state.getOutState(i)))
                strout.append(" with probability " + str(ghmmwrapper.double_array_getitem(state.out_a,i)))
            strout.append( "\nIngoing transitions:")
            for i in range(state.in_states):
                strout.append( "\ntransition from state " + str( state.getInState(i)))
                strout.append( " with probability " + str(ghmmwrapper.double_array_getitem(state.in_a,i)))
            strout.append( "\nint fix:" + str(state.fix) + "\n")

        if self.hasFlags(kSilentStates):
            strout.append("\nSilent states: \n")
            for k in range(hmm.N):
                strout.append( str(self.cmodel.getSilent(k)) + ", ")
        strout.append( "\n")
        return join(strout,'')



    def extendDurations(self, durationlist):
        """ extend states with durations larger one
            this done by explicit state copying in C """

        for i in range(len(durationlist)):
            if durationlist[i] > 1:
                error = self.cmodel.duration_apply(i, durationlist[i])
                if error:
                    log.error( "durations not applied")
                else:
                    self.N = self.cmodel.N

    def getEmission(self, i):
        state = self.cmodel.getState(i)
        if self.hasFlags(kHigherOrderEmissions):
            order = ghmmwrapper.int_array_getitem(self.cmodel.order, i)
            emissions = ghmmhelper.double_array2list(state.b, self.M**(order+1))
        else:                
            emissions = ghmmhelper.double_array2list(state.b, self.M)
        return emissions

    def setEmission(self, i, distributionParameters):
        """ Set the emission distribution parameters for a discrete model."""
        assert len(distributionParameters) == self.M
        # ensure proper indices XXX InsertError
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = self.cmodel.getState(i)

        # updating silent flag and/or model type if necessary
        if self.hasFlags(kSilentStates):
            if sum(distributionParameters) == 0.0:
                self.cmodel.setSilent(i, 1)
            else:
                self.cmodel.setSilent(i, 0)
                #change model_type and free array if no silent state is left
                if 0 == sum(ghmmhelper.int_array2list(self.cmodel.silent,self.N)):
                    self.clearFlags(kSilentStates)
                    ghmmwrapper.free(self.cmodel.silent)
                    self.cmodel.silent = None
        #if the state becomes the first silent state allocate memory and set the silen flag
        elif sum(distributionParameters) == 0.0:
            self.setFlags(kSilentStates)
            slist = [0]*self.N
            slist[i] = 1
            self.cmodel.silent = ghmmhelper.list2int_array(slist)

        #set the emission probabilities
        ghmmwrapper.free(state.b)
        state.b = ghmmhelper.list2double_array(distributionParameters)


    # XXX Change name?
    def backwardTermination(self, emissionSequence, pybeta, scalingVector):
        """

            Result: the backward log probability of emissionSequence
        """
        seq = emissionSequence.cseq.getSequence(0)
        
        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2double_array(scalingVector)
        
        # alllocating beta matrix
        t = len(emissionSequence)
        cbeta = ghmmhelper.list2double_matrix(pybeta)
        #print cbeta[0]

        # allocating double * for log probability
        log_p = ghmmwrapper.double_array_alloc(1)
        
        error = self.cmodel.backward_termination(seq, t, cbeta[0], cscale, log_p)
        if error == -1:
            log.error("backward finished with -1: EmissionSequence cannot be build.")

        logp = ghmmwrapper.double_array_getitem(log_p, 0)

        # deallocation
        ghmmwrapper.free(log_p)
        ghmmwrapper.free(cscale)
        ghmmwrapper.double_matrix_free(cbeta[0],t)
        return logp


    # XXX Implement in C
    def logprob(self, emissionSequence, stateSequence):
        """ log P[ emissionSequence, stateSequence| m] """
        
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        state = self.cmodel.getState(stateSequence[0])
        emissionProb = ghmmwrapper.double_array_getitem(state.b, emissionSequence[0])
        if emissionProb == 0:
            silent = self.hasFlags(kSilentStates) and self.cmodel.silent[stateSequence[0]]
            if silent == 1:
                emissionProb = 1
            else:
                raise SequenceCannotBeBuild, "first symbol " + str(emissionSequence[i+1]) + " not emitted by state " + str(stateSequence[0])
                        
        logP = math.log(state.pi * emissionProb )
        
        symbolIndex = 1

        try:
            for i in range(len(emissionSequence)-1):
                cur_state = self.cmodel.getState(stateSequence[i])
                next_state = self.cmodel.getState(stateSequence[i+1])
                for j in range(cur_state.out_states):
                    out_id = ghmmwrapper.int_array_getitem(cur_state.out_id, j)
                    if out_id == stateSequence[i+1]:
                        emissionProb = ghmmwrapper.double_array_getitem(next_state.b, emissionSequence[symbolIndex])
                        # print "b["+str(emissionSequence[symbolIndex])+"] in state " + str(stateSequence[i+1]) + " = ",emissionProb
                        symbolIndex += 1
                        if emissionProb == 0:
                            silent = self.hasFlags(kSilentStates) and self.cmodel.getSilent(stateSequence[i+1])
                            if silent == 1:
                                emissionProb = 1
                                symbolIndex -= 1
                            else:
                                raise SequenceCannotBeBuild, "symbol " + str(emissionSequence[i+1]) + " not emitted by state "+ str(stateSequence[i+1])

                        logP += math.log( ghmmwrapper.double_array_getitem(state.out_a,j) * emissionProb)
                        break
        except IndexError:
            pass #XXX Funny
        return logP     

    # XXX Make C defaults available to ghmm.py, use baum_welch_nstep only
    def baumWelch(self, trainingSequences, nrSteps = None, loglikelihoodCutoff = None):
        """ Reestimates the model with the sequence in 'trainingSequences'.
           
            Note that training for models including silent states is not yet supported.

            nrSteps is the maximal number of BW-steps
            loglikelihoodCutoff is the least relative improvement in likelihood with respect to the last iteration 
            required to continue.
           
        """
        if not isinstance(trainingSequences,EmissionSequence) and not isinstance(trainingSequences,SequenceSet):
            raise TypeError, "EmissionSequence or SequenceSet required, got " + str(trainingSequences.__class__.__name__)

        if self.hasFlags(kSilentStates):
            raise NotImplementedError("Sorry, training of models containing silent states not yet supported.")
        else:
            if nrSteps == None:
                self.cmodel.baum_welch(trainingSequences.cseq)
            else:
                assert loglikelihoodCutoff != None
                self.cmodel.baum_welch_nstep(trainingSequences.cseq, nrSteps, loglikelihoodCutoff)

    # XXX Go over names (singular/plural)

    def applyBackground(self, backgroundWeight):
        """Apply the background distribution to the emission probabilities of states which
           have been assigned one (usually in the editor and coded in the XML).
           applyBackground computes a convex combination of the emission probability and
           the background, where the backgroundWeight parameter (within [0,1]) controls
           the background's contribution for each state.
        """
        assert len(backgroundWeight) == self.N, "Argument 'backgroundWeight' does not match number of states."
        
        cweights = ghmmhelper.list2double_array(backgroundWeight)
        result = self.cmodel.background_apply(cweights)
        
        ghmmwrapper.free(cweights)
        if result:
            log.error("applyBackground failed.")
						
    
    def setBackground(self, backgroundObject, stateBackground):
        """ Configure model to use the background distributions in 'backgroundObject'. 
            'stateBackground' is a list of indixes (one for each state) refering to distributions
            in 'backgroundObject'.
            
            Note: values in backgroundObject are deepcopied into model
        """
        
        if not isinstance(backgroundObject,BackgroundDistribution):
            raise TypeError, "BackgroundDistribution required, got " + str(emissionSequences.__class__.__name__)        

        if not type(stateBackground) == list:
            raise TypeError, "list required got "+ str(type(stateBackground))
            
        assert len(stateBackground) == self.N, "Argument 'stateBackground' does not match number of states."

        if self.background != None:
            del(self.background)
            ghmmwrapper.free(self.cmodel.background_id)
        self.cmodel.bp = backgroundObject.cbackground
        self.background = backgroundObject
        self.cmodel.background_id = ghmmhelper.list2int_array(stateBackground)

        # updating model type
        self.setFlags(kBackgroundDistributions)
    
    # XXX Unify next two methods
    def assignAllBackgrounds(self,stateBackground):
        """ Change all the assignments of background distributions to states.
        """
        if not type(stateBackground) == list:
           raise TypeError, "list required got "+ str(type(stateBackground))
        
        assert self.cmodel.background_id is not None, "Error: No backgrounds defined in model."   
        assert len(stateBackground) == self.N, "Error: Number of weigths does not match number of states."
        # check for valid background id
        for d in stateBackground:
            assert d in range(self.cbackground.n), "Error: Invalid background distribution id."  

        for i in range(self.N):
            self.cmodel.background_id[i] = stateBackground[i]

    def assignStateBackground(self, state, backgroundID):
        assert self.cmodel.background_id is not None, "Error: No backgrounds defined in model."   
        # XXX CHeck
        if self.labelDomain.isAdmissable(backgroundID):
            self.cmodel.background_id[state] = backgroundID
        else:
            log.error( str(backgroundID) + " is not contained in labelDomain."  )


    def getBackgroundAssignments(self):
        return ghmmhelper.int_array2list(self.cmodel.background_id, self.N)


    def updateTieGroups(self):
        """ XXX Name uninformative. Average emission probabilities of tied states. """
        assert self.hasFlags(kTiedEmissions) and self.cmodel.tied_to is not None, "cmodel.tied_to is undefined."
        self.cmodel.update_tie_groups()


    def setTieGroups(self, tieList):
        assert len(tieList) == self.N, "Number of entries in tieList is different from number of states."
        
        if self.cmodel.tied_to is None:
            log.debug( "allocating tied_to")
            self.cmodel.tied_to = ghmmhelper.list2int_array(tieList)
            self.setFlags(kTiedEmissions)
        else:
            log.debug( "tied_to already initialized")
            for i in range(self.N):
                self.cmodel.tied_to[i] = tieList[i]


    def removeTiegroups(self):
        ghmmwrapper.free(self.cmodel.tied_to)
        self.cmodel.tied_to = None
        self.clearFlags(kTiedEmissions)
    

    def getTieGroups(self):
        assert self.hasFlags(kTiedEmissions) and self.cmodel.tied_to is not None, "cmodel.tied_to is undefined."
        
        return ghmmhelper.int_array2list(self.cmodel.tied_to, self.N)
    

    def getSilentFlag(self,state):
        if self.hasFlags(kSilentStates):
            return self.cmodel.getSilent(state)
        else:
            return 0


    def normalize(self):
        """ Normalize transition probs, emission probs (if applicable) """

        log.debug( "Normalizing now.")

        i_error = self.cmodel.normalize()
        if i_error == -1:
            log.error("normalization failed")


    def toMatrices(self):
        "Return the parameters in matrix form."
        A = []
        B = []
        pi = []
        if self.hasFlags(kHigherOrderEmissions):
            order = ghmmhelper.int_array2list(self.cmodel.order, self.N)
        else:
            order = [0]*self.N

        for i in range(self.cmodel.N):
            A.append([0.0] * self.N)
            state = self.cmodel.getState(i)
            pi.append(state.pi)
            B.append(ghmmhelper.double_array2list(state.b,self.M ** (order[i]+1)))
            for j in range(state.out_states):
                state_index = ghmmwrapper.int_array_getitem(state.out_id, j)
                A[i][state_index] = ghmmwrapper.double_array_getitem(state.out_a,j)

        return [A,B,pi]


    def isSilent(self,state):
        """ Returns True if 'state' is silent, False otherwise
        
        """
        assert 0 <= state <= self.N-1, "Invalid state index"
        
        if self.hasFlags(kSilentStates) and self.cmodel.silent[state]:
            return True
        else:
            return False    


    def pathPosterior(self, sequence, path):
        """ Returns the log posterior probability for 'path' having generated 'sequence'.
            
            CAVEAT: statePosterior needs to calculate the complete forward and backward matrices. If 
            you are interested in multiple paths it would be more efficient to use the 'posterior' function
            directly and not multiple calls to pathPosterior
        """
        # XXX for silent states things are more complicated -> to be done
        if self.hasFlags(kSilentStates):
            raise NotImplementedError, "Models with silent states not yet supported."

        # checking path validity (XXX too inefficient ?)
        for p in path:
            assert 0 <= p <= self.N-1, "Invalid state index "+str(p)+". Model and path are incompatible"

        # calculate complete posterior matrix
        post = self.posterior(sequence)
        path_posterior = []
        
        if not self.hasFlags(kSilentStates):
            # if there are no silent states things are straightforward
            
            assert len(path) == len(sequence), "Path and sequence have different lengths"
            
        
            # appending posteriors for each element of path
            for p in range(len(path)):
                #print "post[",p,"],[",path[p],"] =",post[p][path[p]]
            
                path_posterior.append(post[p][path[p]])
               
            return path_posterior  

        # XXX silent states are yet to be done
#        else:
#            # for silent state models we have to propagate the silent states in each column of the 
#            # posterior matrix 
#            
#            assert not self.isSilent(path[0]), "First state in path must not be silent."
#            
#            j = 0   # path index
#            for i in range(len(sequence)):
#                pp = post[i][path[j]]
#                
#                print pp
#                
#                if pp == 0:
#                    return float('-inf')
#                else:
#                    path_log_lik += math.log(post[p][path[p]])
#                    j+=1
#                
#                
#                # propagate path up until the next emitting state
#                while self.isSilent(path[j]):
#                    
#                    print "** silent state ",path[j]
#                    
#                    pp =  post[i][path[j]]
#                    if pp == 0:
#                        return float('-inf')
#                    else:
#                        path_log_lik += math.log(post[p][path[p]])
#                        j+=1    
#                
#            return path_log_lik        


    def statePosterior(self, sequence, state, time):
        """ Return the log posterior probability for being at 'state' at time 'time' in 'sequence'.
            
            CAVEAT: statePosterior needs to calculate the complete forward and backward matrices. If 
            you are interested in multiple states it would be more efficient to use the hmm.posterior() method
            directly and not multiple calls to statePosterior
        
        """
        # XXX for silent states things arr more complicated -> to be done
        if self.hasFlags(kSilentStates):
            raise NotImplementedError, "Models with silent states not yet supported."

        # checking function arguments
        assert isinstance(sequence, EmissionSequence), "Input to posterior must be EmissionSequence object"
        assert 0 <= time <= len(sequence), "Invalid sequence index: "+str(time)+" (sequence has length "+str(len(sequence))+" )."
        assert 0 <= state <= self.N-1, "Invalid state index: " +str(state)+ " (models has "+str(self.N)+" states )."

        post = self.posterior(sequence)
        return post[time][state]


    def posterior(self, sequence):
        """ Posterior distribution matrix for 'sequence'.
        
        """

        # XXX for silent states things are more complicated -> to be done    	
        if self.hasFlags(kSilentStates):
            raise NotImplementedError, "Models with silent states not yet supported."

        assert isinstance(sequence, EmissionSequence), "Input to posterior must be EmissionSequence object"
        
        (alpha,scale)  = self.forward(sequence)
        beta = self.backward(sequence,scale)
       
        post_mat = []
        for i in range(len(sequence)):
            post = []            
            for state in range(self.N):
                #s = sum(alpha[i])
                #print alpha[i]
                post.append( alpha[i][state] * beta[i][state] ) 
       
            post_mat.append(post)   

        return post_mat


    def write(self,fileName):
        """ Writes HMM to file 'fileName'.

        """
        if self.cmodel.alphabet is None:
            self.cmodel.alphabet = self.emissionDomain.toCstruct()

        self.cmodel.write_xml(fileName)


######################################################
class StateLabelHMM(DiscreteEmissionHMM):
    """ Labelled HMMs with discrete emissions. 
        Same feature list as in DiscreteEmission models.    
    
    """
    def __init__(self, emissionDomain, distribution, labelDomain, cmodel):
        DiscreteEmissionHMM.__init__(self, emissionDomain, distribution, cmodel)

        assert isinstance(labelDomain, LabelDomain), "Invalid labelDomain"
        self.labelDomain = labelDomain


    def __str__(self):
        hmm = self.cmodel
        strout = [str(self.__class__.__name__)]
        if self.cmodel.name:
            strout.append( " " + str(self.cmodel.name))
        strout.append(  "(N= "+ str(hmm.N))
        strout.append(  ", M= "+ str(hmm.M)+')\n')

        f = lambda x: "%.2f" % (x,) # float rounding function 
        
        # XXX 16
        if self.cmodel.model_type &  16: #kHigherOrderEmissions
            order = ghmmhelper.int_array2list(self.cmodel.order, self.N)
        else:
            order = [0]*hmm.N
        label = ghmmhelper.int_array2list(hmm.label, self.N)
        
        if hmm.N <= 4:
            iter_list = range(self.N)
        else:
            iter_list = [0,1,'X',hmm.N-2,hmm.N-1]            

        for k in iter_list:
            if k == 'X':
                strout.append('\n  ...\n\n')
                continue
        
            state = hmm.getState(k)
            strout.append( "  state "+ str(k) +' (')
            if order[k] > 0:
                strout.append( 'order= '+ str(order[k])+',')

            strout.append( "initial= " + f(state.pi)+', label= ' + str(self.labelDomain.external(label[k])) + ')\n')
            strout.append( "    Emissions: ")
            for outp in range(hmm.M**(order[k]+1)):
                strout.append(f(ghmmwrapper.double_array_getitem(state.b,outp)))
                if outp < hmm.M**(order[k]+1)-1:
                    strout.append( ', ')
                else:    
                    strout.append('\n')

            strout.append( "    Transitions:")
            #trans = [0.0] * hmm.N
            for i in range( state.out_states):
                strout.append( " ->" + str( state.getOutState(i))+' ('+ f(ghmmwrapper.double_array_getitem(state.out_a,i) ) +')' )
                if i < state.out_states-1:
                    strout.append( ',')
                #strout.append(" with probability " + str(ghmmwrapper.double_array_getitem(state.out_a,i)))

            strout.append('\n')

        return join(strout,'')


    def verboseStr(self):
        hmm = self.cmodel
        strout = ["\nGHMM Model\n"]
        strout.append("Name: " + str(self.cmodel.name))
        strout.append("\nModelflags: "+ self.printtypes(self.cmodel.model_type))
        strout.append("\nNumber of states: "+ str(hmm.N))
        strout.append("\nSize of Alphabet: "+ str(hmm.M))

        if hmm.model_type & kHigherOrderEmissions:
            order = ghmmhelper.int_array2list(hmm.order, self.N)
        else:
            order = [0]*hmm.N
        label = ghmmhelper.int_array2list(hmm.label, self.N)
        for k in range(hmm.N):
            state = hmm.getState(k)
            strout.append("\n\nState number "+ str(k) +":")

            strout.append("\nState label: "+str(self.labelDomain.external(label[k])))

            strout.append("\nState order: " + str(order[k]))
            strout.append("\nInitial probability: " + str(state.pi))
            strout.append("\nOutput probabilites:\n")
            for outp in range(hmm.M**(order[k]+1)):
                strout+=str(ghmmwrapper.double_array_getitem(state.b,outp))
                if outp % hmm.M == hmm.M-1:
                    strout.append("\n")
                else:
                    strout.append(", ")

            strout.append("Outgoing transitions:")
            for i in range( state.out_states):
                strout.append("\ntransition to state " + str(state.getOutState(i)) + " with probability " + str(state.getOutProb(i)))
            strout.append( "\nIngoing transitions:")
            for i in range(state.in_states):
                strout.append( "\ntransition from state " + str(state.getInState(i)) + " with probability " + str(state.getInProb(i)))
            strout.append("\nint fix:" + str(state.fix) + "\n")

        if hmm.model_type & kSilentStates:
            strout.append("\nSilent states: \n")
            for k in range(hmm.N):
                strout.append(str(hmm.silent[k]) + ", ")
            strout.append("\n")

        return join(strout,'')

    def setLabels(self, labelList):
        """  Set the state labels to the values given in labelList.
             LabelList is in external representation.
        """
        
        assert len(labelList) == self.N, "Invalid number of labels."
        
        # set state label to to the appropiate index
        for i in range(self.N):
            if not self.labelDomain.isAdmissable(labelList[i]):
                raise GHMMOutOfDomain, "Label "+str(labelList[i])+" not included in labelDomain."
            
        ghmmwrapper.free(self.cmodel.label)
        self.cmodel.label = ghmmhelper.list2int_array([self.labelDomain.internal(l) for l in labelList])

    def getLabels(self):
        labels = ghmmhelper.int_array2list(self.cmodel.label, self.N)
        return [self.labelDomain.external(l) for l in labels]
    
    def getLabel(self,stateIndex):
        """ Returns label of the state 'stateIndex'.
        
        """
        return self.cmodel.getLabel(stateIndex)
         
    def externalLabel(self, internal):
        """ Returns label representation of an int or list of int
        """

        if type(internal) is int:
            return self.labelDomain.external[internal] # return Label
        elif type(internal) is list:
            return self.labelDomain.externalSequence(internal)
        else:
            raise TypeError, 'int or list needed'

    def internalLabel(self, external):
        """ Return int representation of an label or list of labels
        """

        if type(external) is list:
            return self.labelDomain.internalSequence(external)
        else:
            return self.labelDomain.internal(external)

    def sampleSingle(self, seqLength, seed = 0):
        seqPtr = self.cmodel.label_generate_sequences(seed, seqLength, 1, seqLength)
        return EmissionSequence(self.emissionDomain, seqPtr, labelDomain = self.labelDomain )

    def sample(self, seqNr,seqLength, seed = 0):
        seqPtr = self.cmodel.label_generate_sequences(seed, seqLength, seqNr, seqLength)
        return SequenceSet(self.emissionDomain,seqPtr, labelDomain = self.labelDomain)


    def viterbiLabels(self, emissionSequences):
        """ 
        Returns the most likely labeling of the input sequence(s) as given by the viterbi path.
        """

        emissionSequences = emissionSequences.asSequenceSet()
        seqNumber = len(emissionSequences)

        assert emissionSequences.emissionDomain == self.emissionDomain, "Sequence and model emissionDomains are incompatible."
        
        (vPath, log_p) = self.viterbi(emissionSequences)

        labels = []
        f = lambda i: self.labelDomain.external(self.getLabel(i))
        if seqNumber == 1:
            labels = map(f, vPath)
        else:
            for j in range(seqNumber):
                labels.append( map(f, vPath[j]) )

        return (labels, log_p)


    def kbest(self, emissionSequences, k = 1):
        """ Compute the k probable labeling for each sequence in emissionSequences

            emissionSequences can either be a SequenceSet or an EmissionSequence

            Result: [l_0, ..., l_T] the labeling of emissionSequences is an emmissionSequence
            object, [[l_0^0, ..., l_T^0], ..., [l_0^k, ..., l_T^k]} for a k-sequence
                    SequenceSet
        """
        if self.hasFlags(kSilentStates):
            raise NotimplementedError("Sorry, k-best decoding on models containing silent states not yet supported.")
        else:
            emissionSequences = emissionSequences.asSequenceSet()
            seqNumber = len(emissionSequences)

            log_p = ghmmwrapper.double_array_alloc(1)

            allLogs = []
            allLabels = []

            for i in range(seqNumber):
                seq = emissionSequences.cseq.getSequence(i)
                seq_len = emissionSequences.cseq.getLength(i)

                labeling = self.cmodel.label_kbest(seq, seq_len, k, log_p)
                oneLabel = ghmmhelper.int_array2list(labeling, seq_len)

                allLabels.append(oneLabel)
                allLogs.append(ghmmwrapper.double_array_getitem(log_p, 0))
                ghmmwrapper.free(labeling)

            ghmmwrapper.free(log_p)

            if emissionSequences.cseq.seq_number > 1:
                return (map(self.externalLabel, allLabels), allLogs)
            else:
                return (self.externalLabel(allLabels[0]), allLogs[0])


    def gradientSearch(self, emissionSequences, eta=.1, steps=20):
        """ trains a model with given sequencesgradescentFunction using gradient descent

            emission_sequences can either be a SequenceSet or an EmissionSequence

        """
        
        # check for labels
        if not self.hasFlags(kLabeledStates):
            raise NotImplementedError("Error: Model is no labeled states.")
        
        emissionSequences = emissionSequences.asSequenceSet()
        seqNumber = len(emissionSequences)

        tmp_model = self.cmodel.label_gradient_descent(emissionSequences.cseq, eta, steps)
        if tmp_model is None:
            log.error("Gradient descent finished not successfully.")
            return False
        else:
            self.cmodel = tmp_model
            return True

    def labelSeqLikelihoods(self, emissionSequences):
        """ Compute a vector ( log( P[s,l| model]) )_{s} of log-likelihoods of the
            individual emission_sequences using the forward algorithm

            emission_sequences is of type SequenceSet

            Result: log( P[emissionSequences,labels| model]) of type float
                    (numarray) vector of floats

        """
        emissionSequences = emissionSequences.asSequenceSet()
        seqNumber = len(emissionSequences)

        assert emissionSequences.cseq.state_labels is not None, "Sequence needs to be labeled."

        likelihood = ghmmwrapper.double_array_alloc(1)
        likelihoodList = []

        for i in range(seqNumber):
            seq = emissionSequences.cseq.getSequence(i)
            labels = ghmmwrapper.get_col_pointer_int(emissionSequences.cseq.state_labels,i)
            tmp = emissionSequences.cseq.getLength(i)
            ret_val = self.cmodel.label_logp(seq, labels, tmp, likelihood)

            if ret_val == -1:
                log.warning("forward returned -1: Sequence"+ str(i) +"cannot be build.")
                likelihoodList.append(-float('Inf'))
            else:
                likelihoodList.append(ghmmwrapper.double_array_getitem(likelihood,0))

        del likelihood
        likelihood = None
        return likelihoodList

    def forwardLabels(self, emissionSequence, labelSequence):
        """

            Result: the (N x T)-matrix containing the forward-variables
                    and the scaling vector
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        logP = ghmmwrapper.double_array_alloc(1)
        n_states = self.cmodel.N

        t = emissionSequence.cseq.getLength(0)
        if t != len(labelSequence):
            raise TypeError, "ERROR: Observation and Labellist must have same length"

        calpha = ghmmwrapper.double_matrix_alloc(t, n_states)
        cscale = ghmmwrapper.double_array_alloc(t)

        seq = emissionSequence.cseq.getSequence(0)
        label = ghmmwrapper.int_array_alloc(t)

        for i in range(len(labelSequence)):
            ghmmwrapper.int_array_setitem(label, i, self.internalLabel(labelSequence[i]))

        error = self.cmodel.label_forward(seq, label, t, calpha, cscale, logP)
        if error == -1:
            log.error( "Forward finished with -1: Sequence " + str(i) + " cannot be build.")

        # translate alpha / scale to python lists
        pyscale = ghmmhelper.double_array2list(cscale, t)
        pyalpha = ghmmhelper.double_matrix2list(calpha,t,n_states)
        logpval = ghmmwrapper.double_array_getitem(logP, 0)
        
        ghmmwrapper.free(label)
        ghmmwrapper.free(logP)
        ghmmwrapper.free(cscale)
        ghmmwrapper.double_matrix_free(calpha,t)
        logP = None
        cscale = None
        calpha = None
        return (logpval, pyalpha, pyscale)

    def backwardLabels(self, emissionSequence, labelSequence, scalingVector):
        """

            Result: the (N x T)-matrix containing the backward-variables
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        t = emissionSequence.cseq.getLength(0)
        if t != len(labelSequence):
            raise TypeError, "ERROR: Observation and Labellist must have same length"

        seq = emissionSequence.cseq.getSequence(0)
        label = ghmmwrapper.int_array_alloc(t)

        for i in range(len(labelSequence)):
            ghmmwrapper.int_array_setitem(label, i, self.internalLabel(labelSequence[i]))

        logP = ghmmwrapper.double_array_alloc(1)

        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2double_array(scalingVector)

        # alllocating beta matrix
        cbeta = ghmmwrapper.double_matrix_alloc(t, self.cmodel.N)

        error = self.cmodel.label_backward(seq, label, t, cbeta, cscale, logP)
        if error == -1:
            log.error( "backward finished with -1: EmissionSequence cannot be build.")

        pybeta = ghmmhelper.double_matrix2list(cbeta,t,self.cmodel.N)
        logpval = ghmmwrapper.double_array_getitem(logP, 0)
        
        # deallocation
        ghmmwrapper.free(logP)
        ghmmwrapper.free(cscale)
        ghmmwrapper.free(label)
        ghmmwrapper.double_matrix_free(cbeta,t)
        cscale = None
        label = None
        cbeta = None
        return (logpval, pybeta)

    def baumWelchLabels(self, trainingSequences, nrSteps = None, loglikelihoodCutoff = None):
        """ Reestimates the model with the sequence in 'trainingSequences'.

            Note that training for models including silent states is not yet supported.

            nrSteps is the maximal number of BW-steps
            loglikelihoodCutoff is the least relative improvement in likelihood with respect to the last iteration
            required to continue.

        """
        if not isinstance(trainingSequences,EmissionSequence) and not isinstance(trainingSequences,SequenceSet):
            raise TypeError, "EmissionSequence or SequenceSet required, got " + str(trainingSequences.__class__.__name__)

        if self.hasFlags(kSilentStates):
            raise NotImplementedError("Sorry, training of models containing silent states not yet supported.")
        else:
            if nrSteps == None:
                self.cmodel.label_baum_welch(trainingSequences.cseq)
            else:
                assert loglikelihoodCutoff != None
                self.cmodel.label_baum_welch_nstep(trainingSequences.cseq, nrSteps, loglikelihoodCutoff)


    def write(self,fileName):
        """ Writes HMM to file 'fileName'.

        """
        if self.cmodel.alphabet is None:
            self.cmodel.alphabet = self.emissionDomain.toCstruct()
            
        if self.cmodel.label_alphabet is None:
            self.cmodel.label_alphabet = self.labelDomain.toCstruct()

        self.cmodel.write_xml(fileName)



class GaussianEmissionHMM(HMM):
    """ HMMs with Gaussian distribution as emissions.
        
    """

    def __init__(self, emissionDomain, distribution, cmodel):
        HMM.__init__(self, emissionDomain, distribution, cmodel)

        # Baum Welch context, call baumWelchSetup to initalize
        self.BWcontext = ""

    def getTransition(self, i, j):
        """ Returns the probability of the transition from state i to state j.
             Raises IndexError if the transition is not allowed
        """
        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."
        assert 0 <= j < self.N, "Index " + str(j) + " out of bounds."
 
     	transition = self.cmodel.get_transition(i, j, 0)
        if transition < 0.0: # Tried to access non-existing edge:
            transition = 0.0
        return transition

    def setTransition(self, i, j, prob):
        """ Accessor function for the transition a_ij """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."
        assert 0 <= j < self.N, "Index " + str(j) + " out of bounds."

        self.cmodel.set_transition(i, j, 0, float(prob))

    def getEmission(self, i):
        """ Return (mu, sigma^2)  """

        # ensure proper index
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = self.cmodel.getState(i)
        mu    = state.getMean(0)
        sigma = state.getStdDev(0)
        return (mu, sigma)

    def setEmission(self, i, (mu, sigma)):
        """ Set the emission distributionParameters for state i """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = self.cmodel.getState(i)
        state.setMean(0, float(mu))  # GHMM C is german: mue instead of mu
        state.setStdDev(0, float(sigma))

    def getMixtureFix(self,state):
        s = self.cmodel.getState(state)
        return ghmmhelper.int_array2list(s.mixture_fix,self.M)
        
        
    def setMixtureFix(self, state ,flags):    
        s = self.cmodel.getState(state)        
        s.mixture_fix = ghmmhelper.list2int_array(flags)
    
    def getStateFix(self,state):
        s = self.cmodel.getState(state)
        return s.fix
        
    def setStateFix(self, state ,flag):    
        s = self.cmodel.getState(state)
        s.fix = flag
    
    def __str__(self):
        hmm = self.cmodel
        strout = [str(self.__class__.__name__)]
        if self.cmodel.name:
            strout.append( " " + str(self.cmodel.name))
        strout.append(  "(N= "+ str(hmm.N)+')\n')

        f = lambda x: "%.2f" % (x,)  # float rounding function      
        
        if hmm.N <= 4:
            iter_list = range(self.N)
        else:
            iter_list = [0,1,'X',hmm.N-2,hmm.N-1]            

        for k in iter_list:
            if k == 'X':
                strout.append('\n  ...\n\n')
                continue
        
            state = hmm.getState(k)
            strout.append("  state "+ str(k) + " (")
            strout.append( "initial= " + f(state.pi) )
            if self.cmodel.cos > 1:
             strout.append(', cos='+ str(self.cmodel.cos))
            strout.append(", mu= " + f(ghmmwrapper.double_array_getitem(state.mue,0))+', ')
            strout.append("sigma= " + f(ghmmwrapper.double_array_getitem(state.u,0)) )
            strout.append(')\n')



            strout.append( "    Transitions: ")
            if self.cmodel.cos > 1:
                strout.append("\n")
            
            for c in range(self.cmodel.cos):
                if self.cmodel.cos > 1:
                    strout.append('      class: ' + str(c)+ ':'  )
                for i in range( state.out_states):
                    strout.append('->' + str(state.getOutState(i)) + ' (' + f(state.getOutProb(i, c))+')' )
                    if i < state.out_states-1:
                        strout.append( ', ')

                strout.append('\n')

        return join(strout,'')


    def verboseStr(self):
        hmm = self.cmodel
        strout = ["\nHMM Overview:"]
        strout.append("\nNumber of states: " + str(hmm.N))
        strout.append("\nNumber of mixture components: " + str(hmm.M))

        for k in range(hmm.N):
            state = hmm.getState(k)
            strout.append("\n\nState number "+ str(k) + ":")
            strout.append("\nInitial probability: " + str(state.pi) + "\n")

            weight = ""
            mue = ""
            u =  ""

            weight += str(ghmmwrapper.double_array_getitem(state.c,0))
            mue += str(ghmmwrapper.double_array_getitem(state.mue,0))
            u += str(ghmmwrapper.double_array_getitem(state.u,0))

            strout.append("  mean: " + str(mue) + "\n")
            strout.append("  variance: " + str(u) + "\n")
            strout.append("  fix: " + str(state.fix) + "\n")
            
            for c in range(self.cmodel.cos):
                strout.append("\n  Class : " + str(c)                )
                strout.append("\n    Outgoing transitions:")
                for i in range( state.out_states):
                    strout.append("\n      transition to state " + str(state.getOutState(i)) + " with probability = " + str(state.getOutProb(i, c)))
                strout.append("\n    Ingoing transitions:")
                for i in range(state.in_states):
                    strout.append("\n      transition from state " + str(state.getInState(i)) +" with probability = "+ str(state.getInProb(i, c)))


        return join(strout,'')



    # different function signatures require overloading of parent class methods    
    def sample(self, seqNr ,T,seed = -1):
        """ Sample emission sequences 
        
        """
        seqPtr = self.cmodel.generate_sequences(seed,T,seqNr,0,-1) 

        return SequenceSet(self.emissionDomain,seqPtr)


    def sampleSingle(self, T,seed= -1):
        """ Sample a single emission sequence of length at most T.
            Returns a Sequence object.
        """
        seqPtr = self.cmodel.generate_sequences(seed,T,1,0,-1) 

        return EmissionSequence(self.emissionDomain,seqPtr)

    def forward(self, emissionSequence):
        """

            Result: the (N x T)-matrix containing the forward-variables
                    and the scaling vector
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        logP = ghmmwrapper.double_array_alloc(1)
        i = self.cmodel.N


        t = emissionSequence.cseq.getLength(0)
        calpha = ghmmwrapper.double_matrix_alloc (t, i)
        cscale = ghmmwrapper.double_array_alloc(t)

        seq = emissionSequence.cseq.getSequence(0)

        error = self.cmodel.forward(seq, t, None, calpha, cscale, logP)
        if error == -1:
            log.error( "Forward finished with -1: Sequence " + str(seq_nr) + " cannot be build.")
        
        # translate alpha / scale to python lists
        pyscale = ghmmhelper.double_array2list(cscale, t)
        pyalpha = ghmmhelper.double_matrix2list(calpha,t,i)

        del logP
        del cscale
        ghmmwrapper.double_matrix_free(calpha,t)
        logP = None
        cscale = None
        calpha = None
        return (pyalpha,pyscale)

    def backward(self, emissionSequence, scalingVector):
        """

            Result: the (N x T)-matrix containing the backward-variables
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        seq = emissionSequence.cseq.getSequence(0)

        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2double_array(scalingVector)

        # alllocating beta matrix
        t = emissionSequence.cseq.getLength(0)
        cbeta = ghmmwrapper.double_matrix_alloc(t, self.cmodel.N)

        error = self.cmodel.backward(seq,t,None,cbeta,cscale)
        if error == -1:
            log.error( "backward finished with -1: EmissionSequence cannot be build.")
            

        pybeta = ghmmhelper.double_matrix2list(cbeta,t,self.cmodel.N)

        # deallocation
        ghmmwrapper.free(cscale)
        ghmmwrapper.double_matrix_free(cbeta,t)
        cscale = None
        cbeta = None
        return pybeta

    def loglikelihoods(self, emissionSequences): 
        """ Compute a vector ( log( P[s| model]) )_{s} of log-likelihoods of the
            individual emission_sequences using the forward algorithm

            emission_sequences is of type SequenceSet

            Result: log( P[emissionSequences| model]) of type float
                    (numarray) vector of floats

        """
        emissionSequences = emissionSequences.asSequenceSet()
        seqNumber = len(emissionSequences)

        if self.cmodel.cos > 1:
            log.debug( "self.cmodel.cos = " + str( self.cmodel.cos) )
            assert self.cmodel.class_change is not None, "Error: class_change not initialized."

        likelihood = ghmmwrapper.double_array_alloc(1)
        likelihoodList = []

        for i in range(seqNumber):
            seq = emissionSequences.cseq.getSequence(i)
            tmp = emissionSequences.cseq.getLength(i)
            
            if self.cmodel.cos > 1:
                self.cmodel.class_change.k = i
            
            ret_val = self.cmodel.logp(seq, tmp, likelihood)
            if ret_val == -1:
                
                log.warning( "forward returned -1: Sequence"+str(i)+"cannot be build.")
                # XXX Eventually this should trickle down to C-level
                # Returning -DBL_MIN instead of infinity is stupid, since the latter allows
                # to continue further computations with that inf, which causes
                # things to blow up later.
                # cmodel.logp() could do without a return value if -Inf is returned
                # What should be the semantics in case of computing the likelihood of
                # a set of sequences
                likelihoodList.append(-float('Inf'))
            else:
                likelihoodList.append(ghmmwrapper.double_array_getitem(likelihood,0))

        ghmmwrapper.free(likelihood)
        likelihood = None

        # resetting class_change->k to default
        if self.cmodel.cos > 1:
           self.cmodel.class_change.k = -1

        return likelihoodList


    
    def viterbi(self, emissionSequences):
        """ Compute the Viterbi-path for each sequence in emissionSequences

            emission_sequences can either be a SequenceSet or an EmissionSequence

            Result: [q_0, ..., q_T] the viterbi-path of emission_sequences is an emmissionSequence
            object, [[q_0^0, ..., q_T^0], ..., [q_0^k, ..., q_T^k]} for a k-sequence
                    SequenceSet
        """
        emissionSequences = emissionSequences.asSequenceSet()
        seqNumber = len(emissionSequences)

        if self.cmodel.cos > 1:
            log.debug( "self.cmodel.cos = "+ str( self.cmodel.cos))
            assert self.cmodel.class_change is not None, "Error: class_change not initialized."


        log_p = ghmmwrapper.double_array_alloc(1)

        allLogs = []
        allPaths = []
        for i in range(seqNumber):
            if self.cmodel.cos > 1:
                # if emissionSequence is a sequenceSet with multiple sequences, 
                # use sequence index as class_change.k
                if emissionSequences.cseq.seq_number > 1:
                    self.cmodel.class_change.k = i

            seq = emissionSequences.cseq.getSequence(i)
            seq_len = emissionSequences.cseq.getLength(i)

            if seq_len != 0:
                viterbiPath = self.cmodel.viterbi(seq,seq_len,log_p)
            else:
                viterbiPath = None

            onePath = []
            for j in range(seq_len):                
                onePath.append( ghmmwrapper.int_array_getitem(viterbiPath,j) )
            

            allPaths.append(onePath)
            allLogs.append(ghmmwrapper.double_array_getitem(log_p, 0))
        
        ghmmwrapper.free(log_p)
        ghmmwrapper.free(viterbiPath)
        viterbiPath = None
        log_p = None

        # resetting class_change->k to default
        if self.cmodel.cos > 1:
           self.cmodel.class_change.k = -1
            
        if emissionSequences.cseq.seq_number > 1:
            return (allPaths, allLogs)
        else:
            return (allPaths[0], allLogs[0])



    def normalize(self):
        """ Normalize transition probs, emission probs (if applicable) """

        for c in range(self.cmodel.cos):
            for i in range(self.N):
                # normalizing transitions
                state = self.cmodel.getState(i)
                out_a_i = ghmmwrapper.get_col_pointer_d(state.out_a,c)
                
                pSum = 0.0
                stateIds = []
                
                for j in range(state.out_states):
                    stateIds.append(state.out_id[j])
                
                    pSum += ghmmwrapper.double_array_getitem(out_a_i,j)
                    
                for j in range(state.out_states):
                    normP = ghmmwrapper.double_array_getitem(out_a_i,j) / pSum
                    out_a_i[j] = normP # updating out probabilities

                    inState = self.cmodel.getState(stateIds[j])
                    in_a =  ghmmwrapper.get_col_pointer_d(inState.in_a,c)
                    for k in range(inState.in_states):
                        inId = inState.in_id[k]
                        if inId == i:
                            in_a[k] = normP # updating in probabilities
                        

    def baumWelch(self, trainingSequences, nrSteps, loglikelihoodCutoff):
        """ Reestimate the model parameters given the training_sequences.
            Perform at most nr_steps until the improvement in likelihood
            is below likelihood_cutoff

            training_sequences can either be a SequenceSet or a Sequence

            Result: Final loglikelihood
        """
        
        if not isinstance(trainingSequences, SequenceSet) and not isinstance(trainingSequences, EmissionSequence):
            raise TypeError, "baumWelch requires a SequenceSet or EmissionSequence object."
        
        assert self.emissionDomain.CDataType == "double", "Continuous sequence needed."
        
        self.baumWelchSetup(trainingSequences, nrSteps)
        ghmmwrapper.ghmm_cmodel_baum_welch(self.BWcontext)        
        likelihood = ghmmwrapper.double_array_getitem(self.BWcontext.logp, 0)
        #(steps_made, loglikelihood_array, scale_array) = self.baumWelchStep(nrSteps,
        #                                                                    loglikelihoodCutoff)
        self.baumWelchDelete()

        return likelihood

    def baumWelchSetup(self, trainingSequences, nrSteps):
        """ Setup necessary temporary variables for Baum-Welch-reestimation.
            Use baumWelchSetup and baumWelchStep if you want more control
            over the training, compute diagnostics or do noise-insertion

            training_sequences can either be a SequenceSet or a Sequence

            Return: a C-array of type ghmm_c_baum_welch_context
        """
        self.BWcontext = ghmmwrapper.ghmm_cmodel_baum_welch_context()
        self.BWcontext.smo  = self.cmodel
        self.BWcontext.sqd  = trainingSequences.cseq    # copy reference to ghmm_cseq
        self.BWcontext.logp = ghmmwrapper.double_array_alloc(1) # place holder for sum of log-likelihood
        self.BWcontext.eps  = 10e-6
        self.BWcontext.max_iter = nrSteps


    def baumWelchStep(self, nrSteps, loglikelihoodCutoff):
        """ Compute one iteration of Baum Welch estimation.
            Use baum_welch_setup and baum_welch_step if you want more control
            over the training, compute diagnostics or do noise-insertion

            training_sequences can either be a SequenceSet or a Sequence
        """
        pass

    def baumWelchDelete(self):
        """ Delete the necessary temporary variables for Baum-Welch-reestimation """

        ghmmwrapper.free(self.BWcontext.logp)
        self.BWcontext.logp = None
        del self.BWcontext
        self.BWcontext = None


    def setPrior(self, prior):
         self.cmodel.prior = prior

         
    def toMatrices(self):
        "Return the parameters in matrix form."
        A = []
        B = []
        pi = []
        for i in range(self.cmodel.N):
            A.append([0.0] * self.N)
            B.append([0.0] * 2)
            state = self.cmodel.getState(i)
            pi.append(state.pi)

            B[i][0] = ghmmwrapper.double_array_getitem(state.mue,0)
            B[i][1] =  ghmmwrapper.double_array_getitem(state.u,0)

            for j in range(state.out_states):
                state_index = ghmmwrapper.int_array_getitem(state.out_id, j)
                A[i][state_index] = ghmmwrapper.double_matrix_getitem(state.out_a,0,j)

        return [A,B,pi]


    def pathPosterior(self, sequence, path):
        """ Returns the log posterior probability for 'path' having generated 'sequence'.
            
            CAVEAT: statePosterior needs to calculate the complete forward and backward matrices. If 
            you are interested in multiple paths it would be more efficient to use the 'posterior' function
            directly and not multiple calls to pathPosterior
        """
        # checking function arguments
        assert isinstance(sequence, EmissionSequence), "Input to posterior must be EmissionSequence object"
        

        # checking path validity (XXX too inefficient ?)
        for p in path:
            assert 0 <= p <= self.N-1, "Invalid state index "+str(p)+". Model and path are incompatible"

        # calculate complete posterior matrix
        post = self.posterior(sequence)
        path_posterior = []
        
        if not self.hasFlags(kSilentStates):
            # if there are no silent states things are straightforward
            assert len(path) == len(sequence), "Path and sequence have different lengths"
        
            # appending posteriors for each element of path
            for p in range(len(path)):
                #print "post[",p,"],[",path[p],"] =",post[p][path[p]]
            
                path_posterior.append(post[p][path[p]])
               
            return path_posterior  

    def statePosterior(self, sequence, state, time):
        """ Return the log posterior probability for being at 'state' at time 'time' in 'sequence'.
            
            CAVEAT: statePosterior needs to calculate the complete forward and backward matrices. If 
            you are interested in multiple states it would be more efficient to use the posterior function
            directly and not multiple calls to statePosterior
        
        """
                    
        # checking function arguments
        assert isinstance(sequence, EmissionSequence), "Input to posterior must be EmissionSequence object"
        assert 0 <= time <= len(sequence), "Invalid sequence index: "+str(time)+" (sequence has length "+str(len(sequence))+" )."
        assert 0 <= state <= self.N-1, "Invalid state index: " +str(state)+ " (models has "+str(self.N)+" states )."

        post = self.posterior(sequence)
        return post[time][state]



    def posterior(self, sequence):
        """ Posterior distribution matrix for 'sequence'.
        
        """
        
        assert isinstance(sequence, EmissionSequence), "Input to posterior must be EmissionSequence object"
        
        (alpha,scale)  = self.forward(sequence)
        beta = self.backward(sequence,scale)
       
        
        post_mat = []
        for i in range(len(sequence)):
            post = []            
            for state in range(self.N):
                #s = sum(alpha[i])
                #print alpha[i]
                post.append( alpha[i][state] * beta[i][state] ) 
       
            post_mat.append(post)   

       

        return post_mat

# XXX - this class will taken over by ContinuousMixtureHMM
class GaussianMixtureHMM(GaussianEmissionHMM):
    """ HMMs with mixtures of Gaussians as emissions.
        Optional features:
            - fixing mixture components in training
        
    
    """
    
    def __init__(self, emissionDomain, distribution, cmodel):
        GaussianEmissionHMM.__init__(self, emissionDomain, distribution, cmodel)

    def getEmission(self, i, comp):
        """ Return (mu, sigma^2, weight) of component 'comp' in state 'i'  """
        state  = self.cmodel.getState(i)
        mu     = state.getMean(comp)
        sigma  = state.getStdDev(comp)
        weigth = state.getWeight(comp)
        return (mu, sigma, weigth)

    def setEmission(self, i, comp,(mu, sigma, weight)):
        """ Set the emission distributionParameters for component 'comp' in state 'i'. """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = self.cmodel.getState(i)
        state.setMean(comp, float(mu))  # GHMM C is german: mue instead of mu
        state.setStdDev(comp, float(sigma))
        state.setWeight(comp, float(weight))

    def getEmissionProbability(self, value, state):
        return self.cmodel.calc_b(state, value)


    def logprob(self, emissionSequence, stateSequence):
        """ log P[ emissionSequence, stateSequence| m] """
        
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)
                        		
        state = self.cmodel.getState(stateSequence[0])
        emissionProb = self.getEmissionProbability(emissionSequence[0],stateSequence[0])
        
    	if (emissionProb == 0): # zero ??? or some small constant?
            raise SequenceCannotBeBuild, "first symbol " + str(emissionSequence[0]) + " not emitted by state " + str(stateSequence[0])
                        
        logP = math.log(state.pi * emissionProb )
        
        #symbolIndex = 1

        try:        
            for i in range(len(emissionSequence)-1):
                cur_state = self.cmodel.getState(stateSequence[i])
                next_state = self.cmodel.getState(stateSequence[i+1])
		
                for j in range(cur_state.out_states):
                    out_id = cur_state.out_id[j]
                    if out_id == stateSequence[i+1]:
                        emissionProb = self.getEmissionProbability(emissionSequence[i+1],out_id)
                        #symbolIndex += 1
                        if emissionProb == 0:
                            raise SequenceCannotBeBuild, "symbol " + str(emissionSequence[i+1]) + " not emitted by state "+ str(stateSequence[i+1])
                        logP += math.log( ghmmwrapper.double_matrix_getitem(cur_state.out_a,0,j) * emissionProb)
                        break
        except IndexError:
            pass
        return logP

    def getPrior(self):
         return self.cmodel.prior

    def setPrior(self, prior):
         self.cmodel.prior = prior  
    
    def __str__(self):
        hmm = self.cmodel
        strout = [str(self.__class__.__name__)]
        if self.cmodel.name:
            strout.append( " " + str(self.cmodel.name))
        strout.append(  "(N= "+ str(hmm.N)+')\n')
        
        f = lambda x: "%.2f" % (x,)  # float rounding function 
        
        if hmm.N <= 4:
            iter_list = range(self.N)
        else:
            iter_list = [0,1,'X',hmm.N-2,hmm.N-1]            

        for k in iter_list:
            if k == 'X':
                strout.append('\n  ...\n\n')
                continue
        
            state = hmm.getState(k)
            strout.append("  state "+ str(k) + " (")
            strout.append( "initial= " + f(state.pi) )
            if self.cmodel.cos > 1:
                strout.append(', cos='+ str(self.cmodel.cos))
            strout.append(')\n')

            weight = ""
            mue = ""
            u =  ""

            for outp in range(hmm.M):
                weight += f(ghmmwrapper.double_array_getitem(state.c,outp))+", "
                mue += f(ghmmwrapper.double_array_getitem(state.mue,outp))+", "
                u += f(ghmmwrapper.double_array_getitem(state.u,outp))+", "

            strout.append( "    Emissions (")
            strout.append("weights= " + str(weight) + ", ")
            strout.append("mu= " + str(mue) + ", ")
            strout.append("sigma= " + str(u) + ")\n")


            strout.append( "    Transitions: ")
            if self.cmodel.cos > 1:
                strout.append("\n")
            
            for c in range(self.cmodel.cos):
                if self.cmodel.cos > 1:
                    strout.append('      class: ' + str(c)+ ':'  )
                for i in range( state.out_states):
                    strout.append('->' + str(state.getOutState(i)) + ' (' + str(state.getOutProb(i, c))+')' )
                    if i < state.out_states-1:
                        strout.append( ', ')

                strout.append('\n')

        return join(strout,'')


    def verboseStr(self):
        "defines string representation"
        hmm = self.cmodel

        strout = ["\nOverview of HMM:"]
        strout.append("\nNumber of states: "+ str(hmm.N))
        strout.append("\nNumber of output distributions per state: "+ str(hmm.M))

        for k in range(hmm.N):
            state = hmm.getState(k)
            strout.append("\n\nState number "+ str(k) +":")
            strout.append("\nInitial probability: " + str(state.pi))
            strout.append("\n"+ str(hmm.M) + " density function(s):\n")

            weight = ""
            mue = ""
            u =  ""

            for outp in range(hmm.M):
                weight += str(ghmmwrapper.double_array_getitem(state.c,outp))+", "
                mue += str(ghmmwrapper.double_array_getitem(state.mue,outp))+", "
                u += str(ghmmwrapper.double_array_getitem(state.u,outp))+", "

            strout.append("  pdf component weights : " + str(weight) + "\n")
            strout.append("  mean vector: " + str(mue) + "\n")
            strout.append("  variance vector: " + str(u) + "\n")

            for c in range(self.cmodel.cos):
                strout.append("\n  Class : " + str(c)                )
                strout.append("\n    Outgoing transitions:")
                for i in range( state.out_states):
                    strout.append("\n      transition to state " + str(state.getOutState(i)) + " with probability = " + str(state.getOutProb(i, c)))
                strout.append("\n    Ingoing transitions:")
                for i in range(state.in_states):
                    strout.append("\n      transition from state " + str(state.getInState(i)) +" with probability = "+ str(state.getInProb(i, c)))

            strout.append("\nint fix:" + str(state.fix) + "\n")
        return join(strout,'')


    def toMatrices(self):
        "Return the parameters in matrix form."
        A = []
        B = []
        pi = []
        for i in range(self.cmodel.N):
            A.append([0.0] * self.N)
            B.append([])
            state = self.cmodel.getState(i)
            pi.append(state.pi)

            B[i].append( ghmmhelper.double_array2list(state.mue,self.cmodel.M) )
            B[i].append( ghmmhelper.double_array2list(state.u,self.cmodel.M) )
            B[i].append( ghmmhelper.double_array2list(state.c,self.cmodel.M) )

            for j in range(state.out_states):
                state_index = ghmmwrapper.int_array_getitem(state.out_id, j)
                A[i][state_index] = ghmmwrapper.double_matrix_getitem(state.out_a,0,j)

        return [A,B,pi]

    def normalize(self):
        """ Normalize transition probs, emission probs (if applicable) """

        for c in range(self.cmodel.cos):
            for i in range(self.N):
                # normalizing transitions
                state = self.cmodel.getState(i)
                out_a_i = ghmmwrapper.get_col_pointer_d(state.out_a,c)
                
                pSum = 0.0
                stateIds = []
                
                for j in range(state.out_states):
                    stateIds.append(state.out_id[j])
                
                    pSum += ghmmwrapper.double_array_getitem(out_a_i,j)
                    
                for j in range(state.out_states):
                    normP = ghmmwrapper.double_array_getitem(out_a_i,j) / pSum
                    out_a_i[j] = normP # updating out probabilities

                    inState = self.cmodel.getState(stateIds[j])
                    in_a =  ghmmwrapper.get_col_pointer_d(inState.in_a,c)
                    for k in range(inState.in_states):
                        inId = inState.in_id[k]
                        if inId == i:
                            in_a[k] = normP # updating in probabilities

                # normalizing mixture weigths
                pSum = 0.0
                for k in range(self.cmodel.M):
                  pSum += ghmmwrapper.double_array_getitem(state.c,k)
               

                for k in range(self.cmodel.M):
                    v = ghmmwrapper.double_array_getitem(state.c,k)
                    state.c[k] =  v / pSum

                      


class ContinuousMixtureHMM(GaussianMixtureHMM):
    """ HMMs with mixtures of any Continuous Distributions as emissions.
        Optional features:
        - fixing mixture components in training
    """

##     NORMAL = 0
##     NORMAL_RIGHT =1
##     NORMAL_LEFT=3
##     UNIFORM=4
     
    def __init__(self, emissionDomain, distribution, cmodel):
        GaussianEmissionHMM.__init__(self, emissionDomain, distribution, cmodel)
 	
 
    def getEmissionProbability(self,value,st):
        return self.cmodel.calc_b(state, value)       
 
    def getEmission(self, i, comp):
        """ Return the paramenters of component 'comp' in state 'i'
        (type, mu, sigma^2, weight) - for a gaussian component
        (type, mu, sigma^2,  min, weight) - for a left trunc gaussian
        (type, mu, sigma^2, max, weight) - for a right trunc gaussian
        (type, max, mix, weight) - for a uniform
        """
        state  = self.cmodel.getState(i)
        mu     = state.getMean(comp) 
        sigma  = state.getStdDev(comp)
        weigth = state.getWeight(comp)
        type   = state.getDensity(comp)
        if ((type == ghmmwrapper.uniform) or (type == ghmmwrapper.normal)):
            return (type, mu, sigma, weigth)
        else:
            a = state.getTruncate(comp)
            return (type, mu, sigma, a, weigth)

    def setEmission(self, i, comp, type, (mu, sigma, weight, a)):
        """ Set the emission distributionParameters for component 'comp' in state 'i'. """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = self.cmodel.getState(i)
        state.setMean(comp, float(mu))  # GHMM C is german: mue instead of mu
        state.setStdDev(comp, float(sigma))
        state.setWeight(comp, float(weight))
        state.setTruncate(comp, float(a))
        state.setDensity(comp, int(type))
	
    def getEmissionProbability(self, value, state):
        return self.cmodel.calc_b(state, value)
         
    def __str__(self):
        "defines string representation"
        hmm = self.cmodel

        strout = "<ContinuousMixtureHMM with "+str(hmm.N)+" states>"

        return join(strout,'')

    def verboseStr(self):
        "defines string representation"
        hmm = self.cmodel

        strout = "\nOverview of HMM:"
        strout.append("\nNumber of states: "+ str(hmm.N))
        strout.append("\nNumber of output distributions per state: "+ str(hmm.M))

        for k in range(hmm.N):
            state = hmm.getState(k)
            strout.append("\n\nState number "+ str(k) +":")
            strout.append("\nInitial probability: " + str(state.pi))
            strout.append("\n"+ str(hmm.M) + " density function(s):\n")

            weight = ""
            mue = ""
            u =  ""
            a = ""
            density = ""
            c = ""

            for outp in range(hmm.M):
                weight += str(ghmmwrapper.double_array_getitem(state.c,outp))+", "
                mue += str(ghmmwrapper.double_array_getitem(state.mue,outp))+", "
                u += str(ghmmwrapper.double_array_getitem(state.u,outp))+", "
                a += str(ghmmwrapper.double_array_getitem(state.a,outp))+", "
                density += str(ghmmwrapper.double_array_getitem(state.a,outp))+","
                c += str(ghmmwrapper.double_array_getitem(state.c,outp))+","

            strout.append("  pdf component weights : " + str(weight) + "\n")
            strout.append("  mean vector: " + str(mue) + "\n")
            strout.append("  variance vector: " + str(u) + "\n")
            strout.append("  a vector: " + str(u) + "\n")
            strout.append("  densities types: " + str(u) + "\n")
            strout.append("  components weitghs: " + str(c) + "\n")

            for c in range(self.cmodel.cos):
                strout.append("\n  Class : " + str(c)                )
                strout.append("\n    Outgoing transitions:")
                for i in range( state.out_states):
                    strout.append("\n      transition to state " + str(state.out_id[i]) + " with probability = " + str(ghmmwrapper.double_matrix_getitem(state.out_a,c,i)))
                strout.append("\n    Ingoing transitions:")
                for i in range(state.in_states):
                    strout.append("\n      transition from state " + str(state.in_id[i]) +" with probability = "+ str(ghmmwrapper.double_matrix_getitem(state.in_a,c,i)))

            strout.append("\nint fix:" + str(state.fix) + "\n")
        return join(strout,'')



    def getPrior(self):
        return self.cmodel.prior
 
    def toMatrices(self):
        """Return the parameters in matrix form.
           It also returns the density type"""
        A = []
        B = []
        pi = []
        d = [] 
        for i in range(self.cmodel.N):
            A.append([0.0] * self.N)
            B.append([])
            state = self.cmodel.getState(i)
            pi.append(state.pi)
 
            B[i].append( ghmmhelper.double_array2list(state.mue,self.cmodel.M) )
            B[i].append( ghmmhelper.double_array2list(state.u,self.cmodel.M) )
            B[i].append( ghmmhelper.double_array2list(state.c,self.cmodel.M) )
            B[i].append( ghmmhelper.double_array2list(state.a,self.cmodel.M) )
             
            d.append([]);

            d[i].append( ghmmhelper.double_array2list(state.density,self.cmodel.M) )
                      
            for j in range(state.out_states):
                state_index = state.out_id[j]
                A[i][state_index] = ghmmwrapper.double_matrix_getitem(state.out_a,0,j)
 
        return [A,B,pi,d]


def HMMDiscriminativeTraining(HMMList, SeqList, nrSteps = 50, gradient = 0):
    """ """
     
    if len(HMMList) != len(SeqList):
        raise TypeError, 'Inputs not equally long'

   
    inplen = len(HMMList)
    if gradient not in [0, 1]:
       raise UnknownInputType, "TrainingType " + gradient + " not supported."
    
    for i in range(inplen):
        if HMMList[i].emissionDomain.CDataType == "double":
            raise TypeError, 'discriminative training is at the moment only implemented on discrete HMMs'
        #initial training with Baum-Welch
        HMMList[i].baumWelch(SeqList[i], 3, 1e-9)

    HMMArray = ghmmwrapper.dmodel_ptr_array_alloc(inplen)
    SeqArray = ghmmwrapper.dseq_ptr_array_alloc(inplen)
    
    for i in range(inplen):
        ghmmwrapper.dmodel_ptr_array_setitem(HMMArray, i, HMMList[i].cmodel)
        ghmmwrapper.dseq_ptr_array_setitem(SeqArray, i, SeqList[i].cseq)

    ghmmwrapper.ghmm_dmodel_label_discriminative(HMMArray, SeqArray, inplen, nrSteps, gradient)

    for i in range(inplen):
        HMMList[i].cmodel = ghmmwrapper.dmodel_ptr_array_getitem(HMMArray, i)
        SeqList[i].cseq   = ghmmwrapper.dseq_ptr_array_getitem(SeqArray, i)
    
    ghmmwrapper.free(HMMArray)
    ghmmwraper.free(SeqArray)

    return HMMDiscriminativePerformance(HMMList, SeqList)


def HMMDiscriminativePerformance(HMMList, SeqList):

    if len(HMMList) != len(SeqList):
        raise TypeRrror, 'Inputs not equally long'

    inplen = len(HMMList)
    
    single = [0.0] * inplen

    HMMArray = ghmmwrapper.dmodel_ptr_array_alloc(inplen)
    SeqArray = ghmmwrapper.dseq_ptr_array_alloc(inplen)

    for i in range(inplen):
        ghmmwrapper.dmodel_ptr_array_setitem(HMMArray, i, HMMList[i].cmodel)
        ghmmwrapper.dseq_ptr_array_setitem(SeqArray, i, SeqList[i].cseq)

    retval = ghmmwrapper.ghmm_dmodel_label_discrim_perf(HMMArray, SeqArray, inplen)
    
    ghmmwrapper.free(HMMArray)
    ghmmwraper.free(SeqArray)

    return retval
        
########## Here comes all the Pair HMM stuff ##########
class DiscretePairDistribution(DiscreteDistribution):
    """
    A DiscreteDistribution over TWO Alphabets: The discrete distribution
    is parameterized by the vector of probabilities.
    To get the index of the vector that corresponds to a pair of characters
    use the getPairIndex method.

    """

    def __init__(self, alphabetX, alphabetY, offsetX, offsetY):
        """
        construct a new DiscretePairDistribution
        @param alphabetX: Alphabet object for sequence X
        @param alphabetY: Alphabet object for sequence Y
        @param offsetX: number of characters the alphabet of sequence X
        consumes at a time
        @param offsetX: number of characters the alphabet of sequence Y
        consumes at a time
        """
        self.alphabetX = alphabetX
        self.alphabetY = alphabetY
        self.offsetX = offsetX
        self.offsetY = offsetY
        self.prob_vector = None
        self.pairIndexFunction = ghmmwrapper.ghmm_dpmodel_pair

    def getPairIndex(self, charX, charY):
        """
        get the index of a pair of two characters in the probability vector
        (if you use the int representation both values must be ints)
        @param charX: character chain or int representation
        @param charY: character chain or int representation
        @return: the index of the pair in the probability vector
        """
        if (not (type(charX) == type(1) and type(charY) == type(1))):
            if (charX == "-"):
                intX = 0 # check this!
            else:
                intX = self.alphabetX.internal(charX)
            if (charY == "-"):
                intY = 0 # check this!
            else:
                intY = self.alphabetY.internal(charY)
        else:
            intX = charX
            intY = charY
        return self.pairIndexFunction(intX, intY,
                                      len(self.alphabetX),
                                      self.offsetX, self.offsetY)

    def setPairProbability(self, charX, charY, probability):
        """
        set the probability of the [air charX and charY to probability
        @param charX: character chain or int representation
        @param charY: character chain or int representation
        @param probability: probability (0<=float<=1)
        """
        self.prob_vector[self.getPairIndex(charX, charY)] = probability

    def getEmptyProbabilityVector(self):
        """
        get an empty probability vector for this distribution (filled with 0.0)
        @return: list of floats
        """
        length = self.pairIndexFunction(len(self.alphabetX) - 1,
                                        len(self.alphabetY) - 1,
                                        len(self.alphabetX),
                                        self.offsetX, self.offsetY) + 1
        return [0.0 for i in range(length)]

    def getCounts(self, sequenceX, sequenceY):
        """
        extract the pair counts for aligned sequences sequenceX and sequenceY
        @param sequenceX: string for sequence X
        @param sequenceY: strinf for sequence Y
        @return: a list of counts
        """
        counts = self.getEmptyProbabilityVector()
        if (self.offsetX != 0 and self.offsetY != 0):
            assert len(sequenceX) / self.offsetX == len(sequenceY) / self.offsetY
            for i in range(len(sequenceX) / self.offsetX):
                charX = sequenceX[i*self.offsetX:(i+1)*self.offsetX]
                charY = sequenceY[i*self.offsetY:(i+1)*self.offsetY]
                counts[self.getPairIndex(charX, charY)] += 1
            return counts
        elif (self.offsetX == 0 and self.offsetY == 0):
            log.error( "Silent states (offsetX==0 and offsetY==0) not supported")
            return counts
        elif (self.offsetX == 0):
            charX = "-"
            for i in range(len(sequenceY) / self.offsetY):
                charY = sequenceY[i*self.offsetY:(i+1)*self.offsetY]
                counts[self.getPairIndex(charX, charY)] += 1
            return counts
        elif (self.offsetY == 0):
            charY = "-"
            for i in range(len(sequenceX) / self.offsetX):
                charX = sequenceX[i*self.offsetX:(i+1)*self.offsetX]
                counts[self.getPairIndex(charX, charY)] += 1
            return counts

class ComplexEmissionSequence(object):
    """
    A complex emission sequence holds the encoded representations of one
    single sequence. The encoding is done by the emissionDomains. It also links
    to the underlying C-structure.

    Note: ComplexEmissionSequence has to be considered imutable for the moment.
    There are no means to manipulate the sequence positions yet.
    """

    def __init__(self, emissionDomains, sequenceInputs, labelDomain = None, labelInput = None):
        """
        @param emissionDomains: a list of EmissionDomain objects corresponding
        to the list of sequenceInputs
        @param sequenceInputs: a list of sequences of the same length (e.g.
        nucleotides and double values) that will be encoded
        by the corresponding EmissionDomain """
        assert len(emissionDomains) == len(sequenceInputs)
        assert len(sequenceInputs) > 0
        self.length = len(sequenceInputs[0])
        for sequenceInput in sequenceInputs:
            assert self.length == len(sequenceInput)
            
        self.discreteDomains = []
        self.discreteInputs = []
        self.continuousDomains = []
        self.continuousInputs = []
        for i in range(len(emissionDomains)):
            if emissionDomains[i].CDataType == "int":
                self.discreteDomains.append(emissionDomains[i])
                self.discreteInputs.append(sequenceInputs[i])
            if emissionDomains[i].CDataType == "double":
                self.continuousDomains.append(emissionDomains[i])
                self.continuousInputs.append(sequenceInputs[i])

        self.cseq = ghmmwrapper.ghmm_dpseq(self.length,
                                           len(self.discreteDomains),
                                           len(self.continuousDomains))
        
        for i in range(len(self.discreteInputs)):
            internalInput = []
            offset = self.discreteDomains[i].getExternalCharacterLength()
            if (offset == None):
                internalInput = self.discreteDomains[i].internalSequence(self.discreteInputs[i])
            else:
                if (type(self.discreteInputs[i]) == type([])):
                    # we have string sequences with equally large characters so
                    # we can join the list representation
                    self.discreteInputs[i] = ("").join(self.discreteInputs[i])
                    
                for j in range(offset - 1):
                    internalInput.append(-1) # put -1 at the start
                for j in range(offset-1, len(self.discreteInputs[i])):
                    internalInput.append(self.discreteDomains[i].internal(
                        self.discreteInputs[i][j-(offset-1):j+1]))
            pointerDiscrete = self.cseq.get_discrete(i)
            for j in range(len(self)):
                ghmmwrapper.int_array_setitem(pointerDiscrete, j, internalInput[j])
            # self.cseq.set_discrete(i, seq)

        for i in range(len(self.continuousInputs)):
            seq = [float(x) for x in self.continuousInputs[i]]
            seq = ghmmhelper.list2double_array(seq)
            pointerContinuous = self.cseq.get_continuous(i)
            for j in range(len(self)):
                ghmmwrapper.double_array_setitem(pointerContinuous, j, self.continuousInputs[i][j])
            # self.cseq.set_continuous(i, seq)

    def __del__(self):
        """
        Deallocation of C sequence struct.
        """
        del self.cseq
        self.cseq = None

    def __len__(self):
        """
        @return: the length of the sequence.
        """
        return self.length

    def getInternalDiscreteSequence(self, index):
        """
        access the underlying C structure and return the internal
        representation of the discrete sequence number 'index'
        @param index: number of the discrete sequence
        @return: a python list of ints
        """
        int_pointer = self.cseq.get_discrete(index)
        internal = ghmmhelper.int_array2list(int_pointer, len(self))
        int_pointer = None
        return internal
    
    def getInternalContinuousSequence(self, index):
        """
        access the underlying C structure and return the internal
        representation of the continuous sequence number 'index'
        @param index: number of the continuous sequence
        @return: a python list of floats
        """
        d_pointer = self.cseq.get_continuous(index)
        internal = ghmmhelper.double_array2list(d_pointer, len(self))
        d_pointer = None
        return internal
    
    def getDiscreteSequence(self, index):
        """
        get the 'index'th discrete sequence as it has been given at the input
        @param index: number of the discrete sequence
        @return: a python sequence
        """
        return self.discreteInputs[index]

    def __getitem__(self, key):
        """
        get a slice of the complex emission sequence
        @param key: either int (makes no big sense) or slice object
        @return: a new ComplexEmissionSequence containing a slice of the
        original
        """
        domains = []
        for domain in self.discreteDomains:
            domains.append(domain)
        for domain in self.continuousDomains:
            domains.append(domain)
        slicedInput = []
        for input in self.discreteInputs:
            slicedInput.append(input[key])
        for input in self.continuousInputs:
            slicedInput.append(input[key])
        return ComplexEmissionSequence(domains, slicedInput)

    def __str__(self):
        """
        string representation. Access the underlying C-structure and return
        the sequence in all it's encodings (can be quite long)
        @return: string representation
        """
        s = ("<ComplexEmissionSequence>")  
        return s

    def verboseStr(self):
        """
        string representation. Access the underlying C-structure and return
        the sequence in all it's encodings (can be quite long)
        @return: string representation
        """
        s = ("ComplexEmissionSequence (len=%i, discrete=%i, continuous=%i)\n"%
             (self.cseq.length, len(self.discreteDomains),
              len(self.continuousDomains)))
        for i in range(len(self.discreteDomains)):
            s += ("").join([str(self.discreteDomains[i].external(x))
                            for x in self.getInternalDiscreteSequence(i)])
            s += "\n"
        for i in range(len(self.continuousDomains)):
            s += (",").join([str(self.continuousDomains[i].external(x))
                            for x in self.getInternalContinuousSequence(i)])
            s += "\n"
        return s
                
class PairHMM(HMM):
    """
    Pair HMMs with discrete emissions over multiple alphabets.
    Optional features: continuous values for transition classes
    """
    def __init__(self, emissionDomains, distribution, cmodel):
        """
        create a new PairHMM object (this should only be done using the
        factory: e.g model = PairHMMOpenXML(modelfile) )
        @param emissionDomains: list of EmissionDomain objects
        @param distribution: (not used) inherited from HMM
        @param cmodel: a swig pointer on the underlying C structure
        """
        HMM.__init__(self, emissionDomains[0], distribution, cmodel)
        self.emissionDomains = emissionDomains
        self.alphabetSizes = []
        for domain in self.emissionDomains:
            if (isinstance(domain, Alphabet)):
                self.alphabetSizes.append(len(domain))

        self.maxSize = 10000
        self.model_type = self.cmodel.model_type  # model type
        self.background = None

        self.states = {}

    def __str__(self):
        """
        string representation (more for debuging) shows the contents of the C
        structure ghmm_dpmodel
        @return: string representation
        """
        hmm = self.cmodel
        strout = ["<PairHMM with "+str(hmm.N)+" states>"]
        return join(strout,'')

    def verboseStr(self):
        """
        string representation (more for debuging) shows the contents of the C
        structure ghmm_dpmodel
        @return: string representation
        """
        hmm = self.cmodel
        strout = ["\nGHMM Model\n"]
        strout.append("Name: " + str(self.cmodel.name))
        strout.append("\nModelflags: "+ self.printtypes(self.cmodel.model_type))
        strout.append("\nNumber of states: "+ str(hmm.N))
        strout.append("\nSize of Alphabet: "+ str(hmm.M))
        for k in range(hmm.N):
            state = hmm.getState(k)
            strout.append("\n\nState number "+ str(k) +":")
            strout.append("\nInitial probability: " + str(state.pi))
            strout.append("\nOutput probabilites: ")
            #strout.append(str(ghmmwrapper.double_array_getitem(state.b,outp)))
            strout.append("\n")

            strout.append("\nOutgoing transitions:")
            for i in range( state.out_states):
                strout.append("\ntransition to state " + str(state.out_id[i]) + " with probability " + str(ghmmwrapper.double_array_getitem(state.out_a,i)))
            strout.append("\nIngoing transitions:")
            for i in range(state.in_states):
                strout.append("\ntransition from state " + str(state.in_id[i]) + " with probability " + str(ghmmwrapper.double_array_getitem(state.in_a,i)))
                strout.append("\nint fix:" + str(state.fix) + "\n")

        if hmm.model_type & kSilentStates:
            strout.append("\nSilent states: \n")
            for k in range(hmm.N):
                strout.append(str(hmm.silent[k]) + ", ")
            strout.append("\n")

        return join(strout,'')


    def viterbi(self, complexEmissionSequenceX, complexEmissionSequenceY):
        """
        run the naive implementation of the Viterbi algorithm and
        return the viterbi path and the log probability of the path
        @param complexEmissionSequenceX: sequence X encoded as ComplexEmissionSequence
        @param complexEmissionSequenceY: sequence Y encoded as ComplexEmissionSequence
        @return: (path, log_p)
        """
        # get a pointer on a double and a int to get return values by reference
        log_p_ptr = ghmmwrapper.double_array_alloc(1)
        length_ptr = ghmmwrapper.int_array_alloc(1)
        # call log_p and length will be passed by reference
        cpath = self.cmodel.viterbi(complexEmissionSequenceX.cseq,
                                    complexEmissionSequenceY.cseq,
                                    log_p_ptr, length_ptr)
        # get the values from the pointers
        log_p = ghmmwrapper.double_array_getitem(log_p_ptr, 0)
        length = length_ptr[0]
        path = [cpath[x] for x in range(length)]
        # free the memory
        ghmmwrapper.free(log_p_ptr)
        ghmmwrapper(length_ptr)
        ghmmwrapper.free(cpath)
        return (path, log_p)
    
    def viterbiPropagate(self, complexEmissionSequenceX, complexEmissionSequenceY, startX=None, startY=None, stopX=None, stopY=None, startState=None, startLogp=None, stopState=None, stopLogp=None):
        """
        run the linear space implementation of the Viterbi algorithm and
        return the viterbi path and the log probability of the path
        @param complexEmissionSequenceX: sequence X encoded as ComplexEmissionSequence
        @param complexEmissionSequenceY: sequence Y encoded as ComplexEmissionSequence
        Optional parameters to run the algorithm only on a segment:
        @param startX: start index in X
        @param startY: start index in Y
        @param stopX: stop index in X
        @param stopY: stop index in Y
        @param startState: start the path in this state
        @param stopState: path ends in this state
        @param startLogp: initialize the start state with this log probability
        @param stopLogp: if known this is the logp of the partial path
        @return: (path, log_p)
        """
        # get a pointer on a double and a int to get return values by reference
        log_p_ptr = ghmmwrapper.double_array_alloc(1)
        length_ptr = ghmmwrapper.int_array_alloc(1)
        # call log_p and length will be passed by reference
        if (not (startX and startY and stopX and stopY and startState and stopState and startLogp)):
            cpath = self.cmodel.viterbi_propagate(
                complexEmissionSequenceX.cseq,
                complexEmissionSequenceY.cseq,
                log_p_ptr, length_ptr,
                self.maxSize)
        else:
            if (stopLogp == None):
                stopLogp = 0
            cpath = self.cmodel.viterbi_propagate_segment(
                complexEmissionSequenceX.cseq,
                complexEmissionSequenceY.cseq,
                log_p_ptr, length_ptr, self.maxSize,
                startX, startY, stopX, stopY, startState, stopState,
                startLogp, stopLogp)
          
        # get the values from the pointers
        log_p = ghmmwrapper.double_array_getitem(log_p_ptr, 0)
        length = length_ptr[0]
        path = [cpath[x] for x in range(length)]
        # free the memory
        ghmmwrapper.free(log_p_ptr)
        ghmmwrapper.free(length_ptr)
        ghmmwrapper.free(cpath)
        return (path, log_p)

    def logP(self, complexEmissionSequenceX, complexEmissionSequenceY, path):
        """
        compute the log probability of two sequences X and Y and a path
        @param complexEmissionSequenceX: sequence X encoded as
        ComplexEmissionSequence
        @param EmissionSequenceEmissionSequenceY: sequence Y encoded as
        ComplexEmissionSequence
        @param path: the state path
        @return: log probability
        """
        cpath = ghmmhelper.list2int_array(path)
        logP = self.cmodel.viterbi_logP(complexEmissionSequenceX.cseq,
                                 complexEmissionSequenceY.cseq,
                                 cpath, len(path))
        ghmmwrapper.free(cpath)
        return logP

    def addEmissionDomains(self, emissionDomains):
        """
        add additional EmissionDomains that are not specified in the XML file.
        This is used to add information for the transition classes.
        @param emissionDomains: a list of EmissionDomain objects
        """
        self.emissionDomains.extend(emissionDomains)
        discreteDomains = []
        continuousDomains = []
        for i in range(len(emissionDomains)):
            if emissionDomains[i].CDataType == "int":
                discreteDomains.append(emissionDomains[i])
                self.alphabetSizes.append(len(emissionDomains[i]))
            if emissionDomains[i].CDataType == "double":
                continuousDomains.append(emissionDomains[i])
                
        self.cmodel.number_of_alphabets += len(discreteDomains)
        self.cmodel.size_of_alphabet = ghmmhelper.list2int_array(self.alphabetSizes)

        self.cmodel.number_of_d_seqs += len(continuousDomains)

    def checkEmissions(self, eps=0.0000000000001):
        """
        checks the sum of emission probabilities in all states
        @param eps: precision (if the sum is > 1 - eps it passes)
        @return: 1 if the emission of all states sum to one, 0 otherwise
        """
        allok = 1
        for state in self.states:
            emissionSum = sum(state.emissions)
            if (abs(1 - emissionSum) > eps):
                log.debug(("Emissions in state %s (%s) do not sum to 1 (%s)" % (state.id, state.label, emissionSum)))
                allok = 0
        return allok

    def checkTransitions(self, eps=0.0000000000001):
        """
        checks the sum of outgoing transition probabilities for all states
        @param eps: precision (if the sum is > 1 - eps it passes)
        @return: 1 if the transitions of all states sum to one, 0 otherwise
        """
        allok = 1
        # from build matrices in xmlutil:
        orders = {}
	k = 0 # C style index
	for s in self.states: # ordering from XML
	    orders[s.index] = k
	    k = k + 1
        
        for state in self.states:
            for tclass in range(state.kclasses):
                outSum = 0.0
                c_state = self.cmodel.getState(orders[state.index])
                for out in range(c_state.out_states):
                    outSum += ghmmwrapper.double_matrix_getitem(c_state.out_a,
                                                        out, tclass)
            
                if (abs(1 - outSum) > eps):
                    log.debug("Outgoing transitions in state %s (%s) do not sum to 1 (%s) for class %s" % (state.id, state.label, outSum, tclass))
                    allok = 0
        return allok
            
class PairHMMOpenFactory(HMMOpenFactory):
    """
    factory to create PairHMM objects from XML files
    """
    def __call__(self, fileName_file_or_dom, modelIndex = None):
        """
        a call to the factory loads a model from a file specified by the
        filename or from a file object or from a XML Document object and
        initializes the model on the C side (libghmm).
        @param fileName_file_or_dom: load the model from a file specified by
        a filename, a file object or a XML Document object
        @param modelIndex: not used (inherited from HMMOpenFactory)
        @return: PairHMM object 
        """
        import xml.dom.minidom
        from ghmm_gato import xmlutil

        if not (isinstance(fileName_file_or_dom, StringIO.StringIO) or
                isinstance(fileName_file_or_dom, xml.dom.minidom.Document)):
            if not os.path.exists(fileName_file_or_dom):
                raise IOError, 'File ' + str(fileName_file_or_dom) + ' not found.'

        hmm_dom = xmlutil.HMM(fileName_file_or_dom)
        if (not hmm_dom.modelType == "pairHMM"):
            raise InvalidModelParameters, "Model type specified in the XML file (%s) is not pairHMM" % hmm_dom.modelType
        # obviously it's a pair HMM
        [alphabets, A, B, pi, state_orders] = hmm_dom.buildMatrices()
        if not len(A) == len(A[0]):
            raise InvalidModelParameters, "A is not quadratic."
        if not len(pi) == len(A):
            raise InvalidModelParameters,  "Length of pi does not match length of A."
        if not len(A) == len(B):
            raise InvalidModelParameters, " Different number of entries in A and B."

        cmodel = ghmmwrapper.ghmm_dp_init()
        cmodel.N = len(A)
        cmodel.M = -1 # no use anymore len(emissionDomain)
        cmodel.prior = -1 # No prior by default
        
        # tie groups are deactivated by default
        cmodel.tied_to = None
        
        # assign model identifier (if specified)
        if hmm_dom.name != None:
            cmodel.name = hmm_dom.name
        else:
            cmodel.name = 'Unused'

        alphabets = hmm_dom.getAlphabets()
        cmodel.number_of_alphabets = len(alphabets.keys())
        sizes = [len(alphabets[k]) for k in alphabets.keys()]
        cmodel.size_of_alphabet = ghmmhelper.list2int_array(sizes)

        # set number of d_seqs to zero. If you want to use them you have to
        # set them manually
        cmodel.number_of_d_seqs = 0

        # c array of states allocated
        cstates = ghmmwrapper.dpstate_array_alloc(cmodel.N)
        # python list of states from xml
        pystates = hmm_dom.state.values()

        silent_flag = 0
        silent_states = []

        maxOffsetX = 0
        maxOffsetY = 0
        
        transitionClassFlag = 0
        maxTransitionIndexDiscrete = len(alphabets.keys())
        maxTransitionIndexContinuous = 0
        
        # from build matrices in xmlutil:
        orders = {}
	k = 0 # C style index
	for s in pystates: # ordering from XML
	    orders[s.index] = k
	    k = k + 1
            
        #initialize states
        for i in range(cmodel.N):
            cstate = ghmmwrapper.dpstate_array_getitem(cstates, i)
            pystate = pystates[i]
            size = len(pystate.itsHMM.hmmAlphabets[pystate.alphabet_id])
            if (pystate.offsetX != 0 and pystate.offsetY != 0):
                size = size**2
            if (len(B[i]) != size):
                raise InvalidModelParameters("in state %s len(emissions) = %i size should be %i" % (pystate.id, len(B[i]), size))
            cstate.b = ghmmhelper.list2double_array(B[i])
            cstate.pi = pi[i]
            if (pi[i] != 0):
                cstate.log_pi = math.log(pi[i])
            else:
                cstate.log_pi = 1

            cstate.alphabet = pystate.alphabet_id
            cstate.offset_x = pystate.offsetX
            cstate.offset_y = pystate.offsetY
            cstate.kclasses = pystate.kclasses
            
            if (pystate.offsetX > maxOffsetX):
                maxOffsetX = pystate.offsetX
            if (pystate.offsetY > maxOffsetY):
                maxOffsetY = pystate.offsetY
                
            if (sum(B[i]) == 0 ): 
                silent_states.append(1)
                silent_flag = 4
            else:
                silent_states.append(0)

                # transition probability
                # cstate.out_states, cstate.out_id, out_a = ghmmhelper.extract_out(A[i])                
                v = pystate.index
                #print "C state index: %i pystate index: %i order: %i" % (i, v, orders[v])
                outprobs = []
                for j in range(len(hmm_dom.G.OutNeighbors(v))):
                    outprobs.append([0.0] * pystate.kclasses)
                myoutid = []
                j = 0                
                for outid in hmm_dom.G.OutNeighbors(v):
                    myorder = orders[outid]
                    myoutid.append(myorder)
                    for tclass in range(pystate.kclasses):
                        outprobs[j][tclass] = hmm_dom.G.edgeWeights[tclass][(v,outid)]
                    j += 1
                cstate.out_states = len(myoutid)
                cstate.out_id = ghmmhelper.list2int_array(myoutid)
                (cstate.out_a, col_len) = ghmmhelper.list2double_matrix(outprobs)
                #set "in" probabilities
                # A_col_i = map( lambda x: x[i], A)
                # Numarray use A[,:i]
                # cstate.in_states, cstate.in_id, cstate.in_a = ghmmhelper.extract_out(A_col_i)
                inprobs = []
                for inid in hmm_dom.G.InNeighbors(v):
                    myorder = orders[inid]
                    # for every class in source
                    inprobs.append([0.0] * pystates[myorder].kclasses)
                myinid = []
                j = 0
                for inid in hmm_dom.G.InNeighbors(v):
                    myorder = orders[inid]
                    myinid.append(myorder)
                    # for every transition class of the source state add a prob
                    for tclass in range(pystates[myorder].kclasses):
                        inprobs[j][tclass] = hmm_dom.G.edgeWeights[tclass][(inid,v)]
                    j += 1
                    
                j = 0
                #for inid in myinid:
                #    print "Transitions (%i, %i)" % (inid ,i)
                #    print inprobs[j]
                #    j += 1
                    
                cstate.in_states = len(myinid)
                cstate.in_id = ghmmhelper.list2int_array(myinid)
                (cstate.in_a, col_len) = ghmmhelper.list2double_matrix(inprobs)
                #fix probabilities by reestimation, else 0
                cstate.fix = 0

                # set the class determination function
                cstate.class_change = ghmmwrapper.ghmm_dp_init_class_change()
                if (pystate.transitionFunction != -1):
                    transitionClassFlag = 1
                    tf = hmm_dom.transitionFunctions[pystate.transitionFunction]
                    # for the moment: do not use the offsets because they
                    # add the risk of segmentation faults at the ends of
                    # the loops or neccessitate index checks at every query
                    # which is not desirable because the transition
                    # functions are used in every iteration. Instead use
                    # shifted input values!
                    if (tf.type == "lt_sum"):
                        ghmmwrapper.set_to_lt_sum(
                            cstate.class_change,
                            int(tf.paramDict["seq_index"]),
                            float(tf.paramDict["threshold"]),
                            0, # int(tf.paramDict["offset_x"]),
                            0) # int(tf.paramDict["offset_y"]))
                        maxTransitionIndexContinuous = max(
                            int(tf.paramDict["seq_index"]),
                            maxTransitionIndexContinuous)
                    elif (tf.type == "gt_sum"):
                        ghmmwrapper.set_to_gt_sum(
                            cstate.class_change,
                            int(tf.paramDict["seq_index"]),
                            float(tf.paramDict["threshold"]),
                            0, # int(tf.paramDict["offset_x"]),
                            0) # int(tf.paramDict["offset_y"]))
                        maxTransitionIndexContinuous = max(
                            int(tf.paramDict["seq_index"]),
                            maxTransitionIndexContinuous)
                    elif (tf.type == "boolean_and"):
                        ghmmwrapper.set_to_boolean_and(
                            cstate.class_change,
                            int(tf.paramDict["seq_index"]),
                            0, # int(tf.paramDict["offset_x"]),
                            0) # int(tf.paramDict["offset_y"]))
                        maxTransitionIndexDiscrete = max(
                            int(tf.paramDict["seq_index"]),
                            maxTransitionIndexDiscrete)
                    elif (tf.type == "boolean_or"):
                        ghmmwrapper.set_to_boolean_or(
                            cstate.class_change,
                            int(tf.paramDict["seq_index"]),
                            0, # int(tf.paramDict["offset_x"]),
                            0) # int(tf.paramDict["offset_y"]))
                        maxTransitionIndexDiscrete = max(
                            int(tf.paramDict["seq_index"]),
                            maxTransitionIndexDiscrete)
                else:
                    ghmmwrapper.ghmm_dp_set_to_default_transition_class(cstate.class_change)

        cmodel.s = cstates

        cmodel.max_offset_x = maxOffsetX
        cmodel.max_offset_y = maxOffsetY
        
        cmodel.model_type += silent_flag
        cmodel.silent = ghmmhelper.list2int_array(silent_states)        
        distribution = DiscreteDistribution(DNA)
        emissionDomains = [Alphabet(hmm_dom.hmmAlphabets[alphabet].name.values()) for alphabet in alphabets]
        model = PairHMM(emissionDomains, distribution, cmodel)
        model.states = pystates
        model.transitionFunctions = hmm_dom.transitionFunctions
        model.usesTransitionClasses = transitionClassFlag
        model.alphabetSizes = sizes
        model.maxTransitionIndexContinuous = maxTransitionIndexContinuous
        model.maxTransitionIndexDiscrete = maxTransitionIndexDiscrete
        return model

PairHMMOpenXML = PairHMMOpenFactory()
