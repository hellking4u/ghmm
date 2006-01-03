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

from Graph import Graph, EdgeWeight
import xmlutil
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

# Initialize logging to stderr
#logging.basicConfig(format="%(asctime)s %(filename)s:%(lineno)d %(levelname)-5s - %(message)s")
log = logging.getLogger("GHMM")

# creating StreamHandler to stderr
hdlr = logging.StreamHandler(sys.stderr)

# setting message format
fmt = logging.Formatter("%(name)s %(asctime)s %(filename)s:%(lineno)d %(levelname)s %(thread)-5s - %(message)s")
hdlr.setFormatter(fmt)

# adding handler to logger object
log.addHandler(hdlr)

# Set the minimal severity of a message to be shown. The levels in
# increasing severity are: DEBUG, INFO, WARNING, ERROR, CRITICAL

log.setLevel(logging.WARNING)
log.info( " I'm the ghmm in "+ __file__)


def logwrapper(level, message):

    if level == 4:
        log.debug(message)
    if level == 3:
        log.info(message)
    if level == 2:
        log.warning(message)
    if level == 1:
        log.error(message)
    if level == 0:
        log.critical(message)


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
        """
        Creates an alphabet out of a listOfCharacters
        @param listOfCharacters: a list of strings (single characters most of
        the time), ints, or other objects that can be used as dictionary keys
        for a mapping of the external sequences to the internal representation
        
        Note: Alphabets should be considered as imutable. That means the
        listOfCharacters and the mapping should never be touched after
        construction.
        """
        self.listOfCharacters = listOfCharacters
        self._lengthOfCharacters = -1
        self.index = {} # Which index belongs to which character
        i = 0
        for c in self.listOfCharacters:
            if (self._lengthOfCharacters != None and type(c) == type("hallo")):
                if (self._lengthOfCharacters == -1):
                    self._lengthOfCharacters = len(c)
                elif (len(c) != self._lengthOfCharacters):
                    self._lengthOfCharacters = None                    
            self.index[c] = i
            i += 1
        if (self._lengthOfCharacters == -1):
            self._lengthOfCharacters = None
        self.CDataType = "int" # flag indicating which C data type should be used

    def __str__(self):
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
        # defining hash and eq is not recommended for mutable types.
        # => listOfCharacters should be considered imutable
        return id(self)

    def size(self):
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
        # this is a gap
        if internal == -1:
            return "-"
        if internal < -1 or len(self.listOfCharacters) < internal:
            raise KeyError, "Internal symbol "+str(internal)+" not recognized."
        return self.listOfCharacters[internal]

    def externalSequence(self, internalSequence):
        """ Given a sequence with the internal representation return the
            external representation
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


DNA = Alphabet(['a','c','g','t'])
AminoAcids = Alphabet(['A','C','D','E','F','G','H','I','K','L',
                       'M','N','P','Q','R','S','T','V','W','Y'])
def IntegerRange(a,b):
    return Alphabet(range(a,b))


# To be used for labelled HMMs. We could use an Alphabet directly but this way it is more explicit.
class LabelDomain(Alphabet):    
    def __init__(self, listOfLabels):
        Alphabet.__init__(self, listOfLabels)


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
class Distribution:
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
        # assert M > index
        # self.mu = mu
        # self.sigma = sigma

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
        assert M > index
        return (self.weigth[i],self.fix[i],self.component[i])

    def check(self):        
        assert self.M == len(self.components)
        assert sum(self.weight) == 1
        assert sum(self.weight > 1) == 0
	assert sum(self.weight < 0) == 0        



#-------------------------------------------------------------------------------
#Sequence, SequenceSet and derived  ------------------------------------------

class EmissionSequence:
    """ An EmissionSequence contains the *internal* representation of
        a sequence of emissions. It also contains a reference to the
        domain where the emission orginated from.
    """

    def __init__(self, emissionDomain, sequenceInput, labelDomain = None, labelInput = None, ParentSequenceSet=None):

        self.emissionDomain = emissionDomain

        if ParentSequenceSet is not None:
            # optional reference to a parent equenceSet. Is needed for reference counting with respect to SequenceSet.__getitem__
            assert isinstance(ParentSequenceSet,SequenceSet), "Error: Invalid reference. Only SequenceSet is valid."    
            self.ParentSequenceSet = ParentSequenceSet
        else:
            self.ParentSequenceSet = None

            
        # check if ghmm is build with asci sequence file support
        if (((isinstance(sequenceInput, str) or isinstance(sequenceInput, unicode)))
            and not ghmmwrapper.ASCI_SEQ_FILE):
            raise UnsupportedFeature ("asci sequence files are deprecated. Please convert your files to the new xml-format or rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        

        if self.emissionDomain.CDataType == "int": # underlying C data type is integer

            # necessary C functions for accessing the ghmm_dseq struct
            self.getPtr = ghmmwrapper.get_col_pointer_int # defines C function to be used to access a single sequence
            self.getSymbol = ghmmwrapper.get_2d_arrayint
            self.setSymbol = ghmmwrapper.set_2d_arrayint
            self.freeFunction = ghmmwrapper.call_ghmm_dseq_free
            self.addSeqFunction = ghmmwrapper.ghmm_dseq_add
            self.freeSubSetFunction = ghmmwrapper.call_ghmm_dseq_subseq_free

            #create a ghmm_dseq with state_labels, if the appropiate parameters are set
            if (isinstance(sequenceInput, list) and (labelInput is not None or labelDomain is not None )):
                assert len(sequenceInput)==len(labelInput)
                assert isinstance(labelInput, list)
                assert isinstance(labelDomain, LabelDomain)
                
                self.labelDomain = labelDomain
                
                internalInput = []
                for i in range(len(sequenceInput)):
                    internalInput.append(self.emissionDomain.internal(sequenceInput[i]))

                #translate the external labels in internal 
                internalLabel = self.labelDomain.internalSequence(labelInput)
                
                (seq,l) = ghmmhelper.list2matrixint([internalInput])
                #make c-array of internal labels
                (label,l) = ghmmhelper.list2matrixint([internalLabel])

                self.cseq = ghmmwrapper.ghmm_dseq_calloc(1)
                #delete allocated space for seq pointers
                ghmmwrapper.freearray(self.cseq.seq)
                self.cseq.seq = seq
                self.cseq.seq_number = 1
                ghmmwrapper.set_arrayint(self.cseq.seq_len,0,l[0])

                #set the state labels and length
                self.cseq.state_labels = label
                self.cseq.state_labels_len = ghmmwrapper.int_array(1)
                ghmmwrapper.set_arrayint(self.cseq.state_labels_len,0,l[0])

            elif isinstance(sequenceInput, list):

                internalInput = map( self.emissionDomain.internal, sequenceInput)
                (seq,l) = ghmmhelper.list2matrixint([internalInput])
                self.cseq = ghmmwrapper.ghmm_dseq_calloc(1)
                #delete allocated space for seq pointers
                ghmmwrapper.freearray(self.cseq.seq)
                self.cseq.seq = seq
                self.cseq.seq_number = 1

                # deactivating labels
                self.cseq.state_labels = None
                self.cseq.state_labels_len = None                
                
                ghmmwrapper.set_arrayint(self.cseq.seq_len,0,l[0]) 

        
            elif (isinstance(sequenceInput, str) or isinstance(sequenceInput, unicode)): # from file

                # reads in the first sequence struct in the input file
                if  not os.path.exists(sequenceInput):
                     raise IOError, 'File ' + str(sequenceInput) + ' not found.'
                else:
                    self.cseq  = ghmmwrapper.seq_read(sequenceInput)

            elif isinstance(sequenceInput, ghmmwrapper.ghmm_dseq):# internal use
                if sequenceInput.seq_number > 1:
                    raise badCPointer, "Use SequenceSet for multiple sequences."
                self.cseq = sequenceInput
                if labelDomain != None:
                    self.labelDomain = labelDomain
                
                
            else:
                raise UnknownInputType, "inputType " + str(type(sequenceInput)) + " not recognized."

        elif self.emissionDomain.CDataType == "double": # underlying C data type is double

            # necessary C functions for accessing the ghmm_cseq struct
            self.getPtr = ghmmwrapper.get_col_pointer_d # defines C function to be used to access a single sequence
            self.getSymbol = ghmmwrapper.get_2d_arrayd
            self.setSymbol = ghmmwrapper.set_2d_arrayd
            self.freeFunction = ghmmwrapper.call_ghmm_cseq_free
            self.addSeqFunction = ghmmwrapper.ghmm_cseq_add
            self.freeSubSetFunction = ghmmwrapper.call_ghmm_cseq_subseq_free
            

            if isinstance(sequenceInput, list):
                (seq,l) = ghmmhelper. list2matrixd([sequenceInput])
                self.cseq = ghmmwrapper.ghmm_cseq_calloc(1)
                ghmmwrapper.freearray(self.cseq.seq)
                self.cseq.seq = seq
                self.cseq.seq_number = 1
                ghmmwrapper.set_arrayint(self.cseq.seq_len,0,l[0])

            elif isinstance(sequenceInput, str) or isinstance(sequenceInput, unicode): # from file
                # reads in the first sequence struct in the input file
                if  not os.path.exists(sequenceInput):
                     raise IOError, 'File ' + str(sequenceInput) + ' not found.'
                else:
                    self.cseq  = ghmmwrapper.seq_d_read(sequenceSetInput)


            elif isinstance(sequenceInput, ghmmwrapper.ghmm_cseq): # internal use
                if sequenceInput.seq_number > 1:
                    raise badCPointer, "Use SequenceSet for multiple sequences."
                self.cseq = sequenceInput
            else:
                raise UnknownInputType, "inputType " + str(type(sequenceInput)) + " not recognized."


        else:
            raise NoValidCDataType, "C data type " + str(self.emissionDomain.CDataType) + " invalid."

        #print "__init__ EmissionSequence: ",self.cseq

    def __del__(self):
        "Deallocation of C sequence struct."
        log.debug( "__del__ EmissionSequence " + str(self.cseq))

        # if a parent SequenceSet exits, we use self.freeSubSetFunction to free memory
        # and clean up the other python-side attributes
        if self.ParentSequenceSet is not None:
            self.ParentSequenceSet = None
            self.emissionDomain = None
            self.freeSubSetFunction(self.cseq)
            # set labelDomain attribute to None, if applicable
            try:
                hasattr(self,labelDomain)
                self.labelDomain = None
            except:
                pass           
            
                
        # otherwise the memory is freed        
        else:
            self.freeFunction(self.cseq)
            self.cseq = None

    def __len__(self):
        "Returns the length of the sequence."
        return ghmmwrapper.get_arrayint(self.cseq.seq_len,0)


    def __setitem__(self, index, value):
        internalValue = self.emissionDomain.internal(value)
        self.setSymbol(self.cseq.seq,0,index,internalValue)
		

    def __getitem__(self, index):
        """ Return the symbol at position 'index'. """
        if index < len(self):
            return self.getSymbol(self.cseq.seq, 0, index)
        else:
            raise IndexError    
    
    def getSeqLabel(self):
        if not ghmmwrapper.SEQ_LABEL_FIELD:
            raise UnsupportedFeature ("the seq_label field is obsolete. If you need it rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        return ghmmwrapper.get_arrayl(self.cseq.seq_label,0)

    def setSeqLabel(self,value):
        if not ghmmwrapper.SEQ_LABEL_FIELD:
            raise UnsupportedFeature ("the seq_label field is obsolete. If you need it rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        ghmmwrapper.set_arrayl(self.cseq.seq_label,0,value)
    
    def getStateLabel(self):
        """ Returns the labeling of the sequence in internal representation"""
        label = []
        if self.cseq.seq_number > index and self.cseq.state_labels != None:
            for j in range(ghmmwrapper.get_arrayint(self.cseq.state_labels_len,0) ):
                label.append(self.labelDomain.external(ghmmwrapper.get_2d_arrayint(self.cseq.state_labels, 0, j)))
            return label
        else:
            raise IndexOutOfBounds(str(0) + " is out of bounds, only " + str(self.cseq.seq_number) + "labels")

    def __str__(self):
        "Defines string representation."
        seq = self.cseq
        strout = []
        strout.append("\nEmissionSequence Instance:\nlength " + str(ghmmwrapper.get_arrayint(seq.seq_len,0))+ ", weight " + str(ghmmwrapper.get_arrayd(seq.seq_w,0))  + ":\n")
        for j in range(ghmmwrapper.get_arrayint(seq.seq_len,0) ):
            strout.append(str( self.emissionDomain.external(self[j]) )   )
            if self.emissionDomain.CDataType == "double":
                strout.append(" ")

        # checking for labels
        if self.emissionDomain.CDataType == "int" and self.cseq.state_labels != None:            
            strout.append("\nState labels:\n")
            for j in range(ghmmwrapper.get_arrayint(seq.state_labels_len,0) ):
                strout.append(str( self.labelDomain.external(ghmmwrapper.get_2d_arrayint(seq.state_labels,0,j)))+ ", ")

    	return join(strout,'')

    def sequenceSet(self):
        """ Return a one-element SequenceSet with this sequence."""
        
        # in order to copy the sequence in 'self', we first create an empty SequenceSet and then
        # add 'self'
        seqSet = SequenceSet(self.emissionDomain, [])
        self.addSeqFunction(seqSet.cseq,self.cseq)
        del(self)
        return seqSet

    def write(self,fileName):
        "Writes the EmissionSequence into file 'fileName'."

        # different function signatures require explicit check for C data type
        if self.emissionDomain.CDataType == "int":
            ghmmwrapper.call_ghmm_dseq_print (fileName, self.cseq)
        if self.emissionDomain.CDataType == "double":
            ghmmwrapper.call_ghmm_cseq_print (fileName, self.cseq,0)

    
    def setWeight(self, value):
        ghmmwrapper.set_arrayd(self.cseq.seq_w,0,value)
        self.cseq.total_w = value
        
    def getWeight(self):
        return ghmmwrapper.get_arrayd(self.cseq.seq_w,0)


class SequenceSet:
    def __init__(self, emissionDomain, sequenceSetInput, labelDomain = None, labelInput = None):
        self.emissionDomain = emissionDomain
        self.cseq = None

        # check if ghmm is build with asci sequence file support
        if ((isinstance(sequenceSetInput, str) or isinstance(sequenceSetInput, unicode))
            and not ghmmwrapper.ASCI_SEQ_FILE):
            raise UnsupportedFeature ("asci sequence files are deprecated. Please convert your files to the new xml-format or rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        
        if self.emissionDomain.CDataType == "int": # underlying C data type is integer
            # necessary C functions for accessing the ghmm_dseq struct
            self.sequenceAllocationFunction = ghmmwrapper.ghmm_dseq_calloc
            self.getPtr = ghmmwrapper.get_col_pointer_int # defines C function to be used to access a single sequence
            self.emptySeq = ghmmwrapper.int_2d_array_nocols # allocate an empty int**
            self.allocSingleSeq = ghmmwrapper.int_array
            self.copySingleSeq = ghmmwrapper.ghmm_dseq_copy
            self.setSeq = ghmmwrapper.set_2d_arrayint_col # assign an int* to a position within an int**
            self.freeFunction = ghmmwrapper.call_ghmm_dseq_free
            self.addSeqFunction = ghmmwrapper.ghmm_dseq_add # add sequences to the underlying C struct
            self.getSymbol = ghmmwrapper.get_2d_arrayint
            self.setSymbolSingle = ghmmwrapper.set_arrayint
            self.getSingleSeq = ghmmwrapper.ghmm_dseq_get_singlesequence
            
            if (isinstance(sequenceSetInput, str)  and labelInput == None): # from file
                # reads in the first sequence struct in the input file
                if  not os.path.exists(sequenceSetInput):
                     raise IOError, 'File ' + str(sequenceSetInput) + ' not found.'
                else:
                     self.cseq  = ghmmwrapper.seq_read(sequenceSetInput)

            #generate a a labeled sequenceSet from a list of lists (sequences), emissionDomain, a list of list (labels) and a model
            elif isinstance(sequenceSetInput, list) and isinstance(labelInput, list) and isinstance(labelDomain, LabelDomain): 
                assert len(sequenceSetInput)==len(labelInput)

                self.labelDomain = labelDomain                
                seq_nr = len(sequenceSetInput)
                self.cseq = ghmmwrapper.ghmm_dseq_calloc(seq_nr)
                self.cseq.seq_number = seq_nr
                
                internalInput = []
                internalLabels = []
                #translate sequences and labels to internal representation
                for i in range(seq_nr):
                    sequenceInput = []
                    for j in range(len(sequenceSetInput[i])):
                        sequenceInput.append(self.emissionDomain.internal(sequenceSetInput[i][j]))

                    internalInput.append(sequenceInput)
                    internalLabels.append( self.labelDomain.internalSequence( labelInput[i]))

                #generate c-arrays
                (seq,lenghts) = ghmmhelper.list2matrixint(internalInput)
                (label,labellen) = ghmmhelper.list2matrixint(internalLabels)

                #delete allocated space for seq pointers
                ghmmwrapper.freearray(self.cseq.seq)
                #set pointers
                self.cseq.seq = seq
                self.cseq.state_labels = label
                self.cseq.state_labels_len = ghmmwrapper.int_array(seq_nr)

                #set lenghts
                for i in range(seq_nr):
                    ghmmwrapper.set_arrayint(self.cseq.seq_len, i, lenghts[i])
                    ghmmwrapper.set_arrayint(self.cseq.state_labels_len, i, labellen[i])

                    
            elif isinstance(sequenceSetInput, list) and labelInput == None:

                seq_nr = len(sequenceSetInput)
                self.cseq = ghmmwrapper.ghmm_dseq_calloc(seq_nr)
                self.cseq.seq_number = seq_nr
                
                internalInput = []
                for i in range(seq_nr):
                    internalInput.append( map( self.emissionDomain.internal, sequenceSetInput[i]))

                (seq,lenghts) = ghmmhelper.list2matrixint(internalInput)
                
                #delete allocated space for seq pointers
                ghmmwrapper.freearray(self.cseq.seq)
                self.cseq.seq = seq
                for i in range(seq_nr):
                    ghmmwrapper.set_arrayint(self.cseq.seq_len,i ,lenghts[i])


            elif isinstance(sequenceSetInput, ghmmwrapper.ghmm_dseq): # inputType == ghmm_dseq*
                self.cseq = sequenceSetInput
                if labelDomain is not None:
                    self.labelDomain = labelDomain
                
            else:    
                raise UnknownInputType, "inputType " + str(type(sequenceSetInput)) + " not recognized."
    
            self.__array = []    
            for index in range(len(self)):
                oneseq = self.getPtr(self.cseq.seq,index)
                self.__array.append(oneseq)


        elif self.emissionDomain.CDataType == "double": # underlying C data type is double
            # necessary C functions for accessing the ghmm_cseq struct
            self.sequenceAllocationFunction =  ghmmwrapper.ghmm_cseq_calloc
            self.getPtr = ghmmwrapper.get_col_pointer_d # defines C function to be used to access a single sequence
            self.emptySeq = ghmmwrapper.double_2d_array_nocols  # cast double* to int** pointer
            self.allocSingleSeq = ghmmwrapper.double_array
            self.copySingleSeq = ghmmwrapper.ghmm_cseq_copy
            self.setSeq = ghmmwrapper.set_2d_arrayd_col # assign a double* to a position within a double**
            self.freeFunction = ghmmwrapper.call_ghmm_cseq_free
            self.addSeqFunction = ghmmwrapper.ghmm_cseq_add # add sequences to the underlying C struct
            self.getSymbol = ghmmwrapper.get_2d_arrayd
            self.setSymbolSingle = ghmmwrapper.set_arrayd
            self.getSingleSeq = ghmmwrapper.ghmm_cseq_get_singlesequence
                        
            if isinstance(sequenceSetInput, list): 
                seq_nr = len(sequenceSetInput)
                self.cseq = ghmmwrapper.ghmm_cseq_calloc(seq_nr)
                self.cseq.seq_number = seq_nr

                (seq,lenghts) = ghmmhelper.list2matrixd(sequenceSetInput)
                #delete allocated space for seq pointers
                ghmmwrapper.freearray(self.cseq.seq)
                self.cseq.seq = seq
                for i in range(seq_nr):
                    ghmmwrapper.set_arrayint(self.cseq.seq_len, i, lenghts[i])

            elif isinstance(sequenceSetInput, str) or isinstance(sequenceSetInput, unicode): # from file
                log.debug( "SequenceSet fromFile" + str (sequenceSetInput))
                if  not os.path.exists(sequenceSetInput):
                     raise IOError, 'File ' + str(sequenceSetInput) + ' not found.'
                else:
                    #self.cseq  = ghmmwrapper.seq_d_read(sequenceSetInput)
                    i = ghmmwrapper.int_array(1)
                    self.__s = ghmmwrapper.ghmm_cseq_read(sequenceSetInput,i)
                    self.cseq = ghmmwrapper.get_seq_d_ptr(self.__s,0)

                                                     
            elif isinstance(sequenceSetInput, ghmmwrapper.ghmm_cseq): # i# inputType == ghmm_cseq**, internal use

                self.cseq = sequenceSetInput
            else:    
                raise UnknownInputType, "inputType " + str(type(sequenceSetInput)) + " not recognized."

            self.__array = []
            for index in range(len(self)):
                oneset = self.getPtr(self.cseq.seq, index)
                self.__array.append(oneset)                

        
        else:
            raise NoValidCDataType, "C data type " + str(self.emissionDomain.CDataType) + " invalid."


    def __del__(self):
        "Deallocation of C sequence struct."
        
        log.debug( "__del__ SequenceSet " + str(self.cseq))
        self.freeFunction(self.cseq)
        self.cseq = None
    
    def __str__(self):
        "Defines string representation."
        seq = self.cseq
        strout =  ["\nNumber of sequences: " + str(seq.seq_number)]
        
        if self.emissionDomain.CDataType == "int":
           for i in range(seq.seq_number):
                strout.append("\nSeq " + str(i)+ ", length " + str(ghmmwrapper.get_arrayint(seq.seq_len,i))+ ", weight " + str(ghmmwrapper.get_arrayd(seq.seq_w,i))  + ":\n")
                for j in range(ghmmwrapper.get_arrayint(seq.seq_len,i) ):
                    strout.append(str( self.emissionDomain.external(( ghmmwrapper.get_2d_arrayint(self.cseq.seq, i, j) )) ))

                # checking for labels 
                if self.emissionDomain.CDataType == "int" and self.cseq.state_labels != None:            
                    strout.append("\nState labels:\n")
                    for j in range(ghmmwrapper.get_arrayint(seq.state_labels_len,i) ):
                        strout.append(str( self.labelDomain.external(ghmmwrapper.get_2d_arrayint(seq.state_labels,i,j))) +", ")

        if self.emissionDomain.CDataType == "double":
            for i in range(seq.seq_number):
                strout.append("\nSeq " + str(i)+ ", length " + str(ghmmwrapper.get_arrayint(seq.seq_len,i))+ ", weight " + str(ghmmwrapper.get_arrayd(seq.seq_w,i))  + ":\n")
                for j in range(ghmmwrapper.get_arrayint(seq.seq_len,i) ):
                    strout.append(str( self.emissionDomain.external(( ghmmwrapper.get_2d_arrayd(self.cseq.seq, i, j) )) ) + " ")

        return join(strout,'')
    


    def __len__(self):
        """ Return the number of sequences in the SequenceSet. """
        return self.cseq.seq_number

    def sequenceLength(self, i):
        """ Return the lenght of sequence 'i' in the SequenceSet """
        return ghmmwrapper.get_arrayint(self.cseq.seq_len,i)
        

    def getWeight(self, i):
        """ Return the weight of sequence i. Weights are used in Baum-Welch"""
        return ghmmwrapper.get_arrayd(self.cseq.seq_w,i)

    def setWeight(self, i, w):
        """ Set the weight of sequence i. Weights are used in Baum-Welch"""
        return ghmmwrapper.set_arrayd(self.cseq.seq_w,i,w)
        
    def __getitem__(self, index):
        """ Return an EmissionSequence object initialized with a reference to 
        sequence 'index'.
        
        """
        # check the index for correct range
        if index >= self.cseq.seq_number:
            raise IndexError
        
        seq = self.getSingleSeq(self.cseq,index) 
        return EmissionSequence(self.emissionDomain, seq, ParentSequenceSet=self) 
    
    
    def getSeqLabel(self,index):
        if not ghmmwrapper.SEQ_LABEL_FIELD:
            raise UnsupportedFeature ("the seq_label field is obsolete. If you need it rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        return ghmmwrapper.get_arrayl(self.cseq.seq_label,index)

    def setSeqLabel(self,index,value):
        if not ghmmwrapper.SEQ_LABEL_FIELD:
            raise UnsupportedFeature ("the seq_label field is obsolete. If you need it rebuild the GHMM with the conditional \"GHMM_OBSOLETE\".")
        ghmmwrapper.set_arrayl(self.cseq.seq_label,index,value)

    def getSequence(self, index):
        """ Returns the index-th sequence in internal representation"""
        seq = []
        if self.cseq.seq_number > index:
            for j in range(ghmmwrapper.get_arrayint(self.cseq.seq_len,index) ):
                seq.append(ghmmwrapper.get_2d_arrayint(self.cseq.seq, index, j))
            return seq
        else:
            raise IndexOutOfBounds(str(index) + " is out of bounds, only " + str(self.cseq.seq_number) + "sequences")

    def getStateLabel(self,index):
        """ Returns the labeling of the index-th sequence in internal representation"""
        label = []
        if self.cseq.seq_number > index and self.cseq.state_labels != None:
            for j in range(ghmmwrapper.get_arrayint(self.cseq.state_labels_len,index) ):
                    label.append(self.labelDomain.external(ghmmwrapper.get_2d_arrayint(self.cseq.state_labels, index, j)))
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
        
        self.addSeqFunction(self.cseq, emissionSequences.cseq)
        del(emissionSequences) # removing merged sequences

    def getSubset(self, seqIndixes):
        """ Returns a SequenceSet containing (references to) the sequences with the indixes in
            'seqIndixes'.
       
        """
        seqNumber = len(seqIndixes)
        seq = self.sequenceAllocationFunction(seqNumber)
        #seq.seq = self.emptySeq(seqNumber)
        
        # checking for state labels in the source C sequence struct
        if self.emissionDomain.CDataType == "int" and self.cseq.state_labels is not None:
            
            log.debug( "SequenceSet: found labels !")
            seq.state_labels = ghmmwrapper.int_2d_array_nocols(seqNumber)
            seq.state_labels_len = ghmmwrapper.int_array(seqNumber)

        for i in range(seqNumber):

            len_i = ghmmwrapper.get_arrayint(self.cseq.seq_len,seqIndixes[i])
            
            #seq_i = self.allocSingleSeq(len_i)
            #source_i = self.getPtr(self.cseq.seq, seqIndixes[i])
            #self.copySingleSeq(seq_i,source_i,len_i)
            
            self.setSeq(seq.seq,i,self.__array[seqIndixes[i]])

            ghmmwrapper.set_arrayint(seq.seq_len,i,len_i)

            # Above doesnt copy seq_id or seq_label or seq_w
            seq_id = int(ghmmwrapper.get_arrayd(self.cseq.seq_id, seqIndixes[i]))
            ghmmwrapper.set_arrayd(seq.seq_id, i, seq_id)
            if ghmmwrapper.SEQ_LABEL_FIELD:
                seq_label = ghmmwrapper.get_arrayl(self.cseq.seq_label, i)
                ghmmwrapper.set_arrayl(seq.seq_label, i, int(seq_label))

            seq_w = ghmmwrapper.get_arrayd(self.cseq.seq_w, i)
            ghmmwrapper.set_arrayd(seq.seq_w, i, seq_w)
             
            # setting labels if appropriate
            # XXX needs to be real copy! XXX
            if self.emissionDomain.CDataType == "int" and self.cseq.state_labels is not None:
                self.setSeq(seq.state_labels, i, ghmmwrapper.get_col_pointer_int(self.cseq.state_labels, seqIndixes[i]))
                ghmmwrapper.set_arrayint(seq.state_labels_len, i, ghmmwrapper.get_arrayint(self.cseq.state_labels_len, seqIndixes[i]))

        seq.seq_number = seqNumber
        
        return SequenceSetSubset(self.emissionDomain, seq, self)
        
    def write(self,fileName):
        "Writes (appends) the SequenceSet into file 'fileName'."
        
        # different function signatures require explicit check for C data type
        if self.emissionDomain.CDataType == "int":
            ghmmwrapper.call_ghmm_dseq_print(fileName, self.cseq)
        if self.emissionDomain.CDataType == "double":    
            ghmmwrapper.call_ghmm_cseq_print(fileName, self.cseq,0)

class SequenceSetSubset(SequenceSet):
    """ 
    SequenceSetSubset contains a subset of the sequences from a SequenceSet object.
    On the C side only the references are used.
    """
    def __init__(self, emissionDomain, sequenceSetInput, ParentSequenceSet , labelDomain = None, labelInput = None):
        # reference on the parent SequenceSet object
        self.ParentSequenceSet =  ParentSequenceSet
        SequenceSet.__init__(self, emissionDomain, sequenceSetInput, labelDomain, labelInput)

        if self.emissionDomain.CDataType == "int":
            self.freeFunction = ghmmwrapper.call_ghmm_dseq_subseq_free
        elif self.emissionDomain.CDataType == "double":
            self.freeFunction = ghmmwrapper.call_ghmm_cseq_subseq_free
    
    def __del__(self):
        """ Since we do not want to deallocate the sequence memory, the destructor has to be
            overloaded.
        
        """
        
        self.freeFunction(self.cseq)
        self.cseq = None
            
        # remove reference on parent SequenceSet object
        self.ParentSequenceSet = None
    


def SequenceSetOpen(emissionDomain, fileName):
    """ Reads a sequence file with multiple sequence sets. 

    Returns a list of SequenceSet objects.
    
    """

    if not os.path.exists(fileName):
        raise IOError, 'File ' + str(fileName) + ' not found.'

    
    if emissionDomain.CDataType == "int":
        readFile = ghmmwrapper.ghmm_dseq_read
        seqPtr = ghmmwrapper.get_seq_ptr
    elif emissionDomain.CDataType == "double":
        readFile = ghmmwrapper.ghmm_cseq_read
        seqPtr = ghmmwrapper.get_seq_d_ptr
            
    dArr = ghmmwrapper.int_array(1)

    structArray = readFile(fileName, dArr)
    setNr = ghmmwrapper.get_arrayint(dArr,0)

    sequenceSets = []
    for i in range(setNr):
        seq = seqPtr(structArray,i)
        sequenceSets.append(SequenceSet(emissionDomain, seq) )

        # setting labels to NULL
        sequenceSets[i].cseq.state_labels = None
        sequenceSets[i].cseq.state_labels_len = None
       
    ghmmwrapper.freearray(dArr)
    return  sequenceSets




#-------------------------------------------------------------------------------
# HMMFactory and derived  -----------------------------------------------------
class HMMFactory:
    """ A HMMFactory is the base class of HMM factories.
        A HMMFactory has just a constructor and a () method
    """


GHMM_FILETYPE_SMO = 'smo'
GHMM_FILETYPE_XML = 'xml'
GHMM_FILETYPE_HMMER = 'hmm'


class HMMOpenFactory(HMMFactory):

    def __init__(self, defaultFileType=None):
        if defaultFileType:
            self.defaultFileType = defaultFileType

    def __call__(self, fileName, modelIndex = None):
        
        if not isinstance(fileName,StringIO.StringIO):
            if not os.path.exists(fileName):
                raise IOError, "File %s not found."% str(fileName)
            
        # XML file: both new and old format
    	if self.defaultFileType == GHMM_FILETYPE_XML:
            try:
                m = self.openNewXML(fileName, modelIndex)
            except WrongFileType:
                m = self.openOldXML(fileName)
                
            return m
	    
        elif self.defaultFileType == GHMM_FILETYPE_SMO:
            return self.openSMO(fileName, modelIndex)

        else:
            # HMMER format models
            return self.openHMMER(fileName)


    def openNewXML(self, fileName, modelIndex):
        # check the type of hmm
        # start the model

        file = ghmmwrapper.parseHMMDocument(fileName)
        if file == None:
            log.debug( "XML has file format problems!")
            raise WrongFileType("file is not in GHMM xml format")

        nrModels = file.noModels
        modelType = file.modelType

        if (modelType == ghmmwrapper.kContinuousHMM):                
            emission_domain = Float()
            distribution = ContinuousMixtureDistribution
            hmmClass = ContinuousMixtureHMM
            getPtr = ghmmwrapper.get_smodel_ptr
            models = ghmmwrapper.get_model_c(file)

        else:
            log.warning("Non-supported model type")

        # for a mixture of HMMs return a list of HMMs
        if nrModels>1 and modelIndex == None:
            result = []
            for i in range(nrModels):
                cmodel = getPtr(models,i)
                m = hmmClass(emission_domain, distribution(emission_domain), cmodel)
                result.append(m)
                
        # for a single HMM just return the HMM
        else:
            if modelIndex == None:
                cmodel = getPtr(models, 0)
            elif modelIndex < nrModels:
                cmodel = getPtr(models, modelIndex)
            else:
                raise IndexOutOfBounds("the file %s has only %s models"% fileName, str(nrModels))

            result = hmmClass(emission_domain, distribution(emission_domain), cmodel)
            

        ghmmwrapper.freearray(models)
        return result

    def openOldXML(self, fileName):
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

            if background_dist != {}:
                 ids = [-1]*m.N
                 for s in hmm_dom.state.values():
                      ids[s.index-1] = s.background # s.index ranges from [1, m.N]  

                 m.setBackground(bg, ids)
                 log.debug( "model_type %x" % m.cmodel.model_type)
                 log.debug("background_id" + str( ghmmhelper.arrayint2list(m.cmodel.background_id, m.N)))
            else:
                 m.cmodel.bp = None
                 m.cmodel.background_id = None

            # check for tied states
            tied = hmm_dom.getTiedStates()
            if len(tied) > 0:
                m.cmodel.model_type += 8  #kTiedEmissions
                m.cmodel.tied_to = ghmmhelper.list2arrayint(tied)

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
        
        nrModelPtr = ghmmwrapper.int_array(1)
    	    
        # XXX broken since silent states are not supported by .smo file format 
        if hmmClass == DiscreteEmissionHMM:
            models = ghmmwrapper.ghmm_d_read(fileName, nrModelPtr)
            getPtr = ghmmwrapper.get_model_ptr
        else:
            models = ghmmwrapper.ghmm_c_read(fileName, nrModelPtr)
            getPtr = ghmmwrapper.get_smodel_ptr

        nrModels = ghmmwrapper.get_arrayint(nrModelPtr, 0)
        if modelIndex == None:
            result = []
            for i in range(nrModels):
                cmodel = getPtr(models, i)
                m = hmmClass(emission_domain, distribution(emission_domain), cmodel)
                result.append(m)
        else:
            if modelIndex < nrModels:
                cmodel = getPtr(models, modelIndex)
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

        [A,B,pi,modelName] = h.getGHMMmatrices()
        return  HMMFromMatrices(emission_domain, distribution, A, B, pi, hmmName=modelName)

        
    def determineHMMClass(self, fileName):
        #
        # smo files 
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
HMMOpen = HMMOpenFactory(GHMM_FILETYPE_SMO)
HMMOpenXML = HMMOpenFactory(GHMM_FILETYPE_XML)


def readMultipleHMMERModels(fileName):
    """
        Reads a file containing multiple HMMs in HMMER format, returns list of
        HMM objects.

    """
    
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
                raise InvalidModelParameters, "Either labelDomain and labelInput are both given or neither of the two."
            
            if isinstance(distribution,DiscreteDistribution):
                
                # HMM has discrete emissions over finite alphabet: DiscreteEmissionHMM
                cmodel = ghmmwrapper.new_model()
                cmodel.N = len(A)
                cmodel.M = len(emissionDomain)
                cmodel.prior = -1 # No prior by default
                
                # tie groups are deactivated by default
                cmodel.tied_to = None
                
                # assign model identifier (if specified)
                if hmmName != None:
                    cmodel.name = hmmName
                else:
                    cmodel.name = 'Unused'

                states = ghmmwrapper.arraystate(cmodel.N)

                silent_flag = 0
                silent_states = []

                tmpOrder = []

                #initialize states
                for i in range(cmodel.N):
                    state = ghmmwrapper.get_stateptr(states,i)
                    # compute state order
                    if cmodel.M > 1:
                        order = math.log(len(B[i]), cmodel.M)-1
                    else:
                        order = len(B[i]) - 1
                        
                    log.debug( "order in state "+str(i)+" = "+ str(order) )
                    # check or valid number of emission parameters
                    order = int(order)
                    if  cmodel.M**(order+1) == len(B[i]):
                        tmpOrder.append(order)
                    else:
                        raise InvalidModelParameters, "The number of "+str(len(B[i]))+ " emission parameters for state "+str(i)+" is invalid. State order can not be determined."
                    
                    state.b = ghmmhelper.list2arrayd(B[i])
                                        
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
                    cmodel.model_type += silent_flag
                    cmodel.silent = ghmmhelper.list2arrayint(silent_states)

                cmodel.maxorder = max(tmpOrder)
                if cmodel.maxorder > 0:
                    log.debug( "Set kHigherOrderEmissions.")
                    cmodel.model_type += 16     #kHigherOrderEmissions
                    cmodel.order = ghmmhelper.list2arrayint(tmpOrder)

                # initialize lookup table for powers of the alphabet size,
                # speeds up models with higher order states
                powLookUp = [1] * (cmodel.maxorder+2)
                for i in range(1,len(powLookUp)):
                    powLookUp[i] = powLookUp[i-1] * cmodel.M
                cmodel.pow_lookup = ghmmhelper.list2arrayint(powLookUp)
                
                # check for state labels
                if labelDomain is not None and labelList is not None:
                    if not isinstance(labelDomain,LabelDomain):
                        raise TypeError, "LabelDomain object required."
                    
                    cmodel.model_type += 64     #kClassLabels
                    m = StateLabelHMM(emissionDomain, distribution, labelDomain, cmodel)
                    m.setLabels(labelList)
                    return m
                else:    
                    return DiscreteEmissionHMM(emissionDomain, distribution, cmodel)

            else:
                raise GHMMError(type(distribution), "Not a valid distribution for Alphabet") 
        else:

            if isinstance(distribution,GaussianDistribution):
                
                cmodel = ghmmwrapper.new_smodel()
                
                # def __call__(self, emissionDomain, distribution, A, B, pi, hmmName = None, labelDomain= None, labelList = None):                


                cmodel.M = 1 # Number of mixture componenent for emission distribution
                cmodel.prior = -1 # Unused

                # determining number of transition classes
                cos = 0
                if type(A[0][0]) == list:
                    cos = len(A)
                    cmodel.N = len(A[0])
                    
                    # allocating class switching context
                    ghmmwrapper.ghmm_c_class_change_alloc(cmodel)
                else: 
                    cos = 1
                    cmodel.N = len(A)
                    A = [A]
                
                cmodel.cos = ghmmhelper.classNumber(A)  # number of transition classes in GHMM

                log.debug( "cmodel.cos = "+str(cmodel.cos))

                states = ghmmwrapper.arraysstate(cmodel.N)

                # Switching functions and transition classes are handled
                # elswhere

                #initialize states
                for i in range(cmodel.N):

                    state = ghmmwrapper.get_sstate_ptr(states,i)
                    state.pi = pi[i]
                    state.M = 1

                    # allocate arrays of emmission parameters
                    state.c = ghmmhelper.list2arrayd([1.0]) # Mixture weights. Unused
                    (mu, sigma) = B[i]
                    state.mue = ghmmhelper.list2arrayd([mu]) #mu = mue in GHMM C-lib.
                    state.u = ghmmhelper.list2arrayd([sigma])
                    state.c = ghmmhelper.list2arrayd([1.0])
                    state.a = ghmmhelper.list2arrayd([0.0])
                    
                    # mixture fixing deactivated by default
                    state.mixture_fix = ghmmhelper.list2arrayint([0])

                    # setting densities types (all normal by default)                    
                    densities = ghmmwrapper.arraydensity(1)
                    state.density = densities
                    ghmmwrapper.set_density(state,0,0)                                
                    
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
                
                cmodel = ghmmwrapper.new_smodel()
                cmodel.M = len(B[0][0]) # Number of mixture componenents for emission distribution
                cmodel.prior = -1 # Unused
                
                # determining number of transition classes
                cos = 0
                if type(A[0][0]) == list:
                    cos = len(A)
                    cmodel.N = len(A[0])
                    # allocating class switching context
                    ghmmwrapper.ghmm_c_class_change_alloc(cmodel)
                    
                else: 
                    cos = 1
                    cmodel.N = len(A)
                    A = [A]
                    
                cmodel.cos = ghmmhelper.classNumber(A)  # number of transition classes in GHMM
                
                states = ghmmwrapper.arraysstate(cmodel.N)

                # Switching functions and transition classes are handled
                # elswhere

                #initialize states
                for i in range(cmodel.N):
                    state = ghmmwrapper.get_sstate_ptr(states,i)
                    state.pi = pi[i]
		    state.M = len(B[0][0])

                    # allocate arrays of emmission parameters
                    state.c = ghmmhelper.list2arrayd([1.0]) # Mixture weights. Unused
                    mu_list = B[i][0]
                    sigma_list = B[i][1]
                    weight_list = B[i][2]
                    
                    state.mue = ghmmhelper.list2arrayd(mu_list) #mu = mue in GHMM C-lib.
                    state.u = ghmmhelper.list2arrayd(sigma_list)
                    state.c = ghmmhelper.list2arrayd(weight_list)
		    state.a = ghmmhelper.list2arrayd([0.0] * state.M)

                    # setting densities types (all normal by default)                    
                    densities = ghmmwrapper.arraydensity(cmodel.M)
                    state.density = densities
                    for j in range(state.M):
                      ghmmwrapper.set_density(state,j,0)

                    # mixture fixing deactivated by default
                    state.mixture_fix = ghmmhelper.list2arrayint([0] * state.M)
                    
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
                 #print " ** mixture model"
                 
                 # Interpretation of B matrix for the mixture case (Example with three states and two components each):
                 #  B = [ 
                 #      [ ["mu11","mu12"],["sig11","sig12"],["a11","a12"],["w11","w12"]   ],
                 #      [  ["mu21","mu22"],["sig21","sig22"],["w21","w22"]  ],
                 #      [  ["mu31","mu32"],["sig31","sig32"],["w31","w32"]  ],
                 #      ]

                assert densities != None, "Continuous Mixture Distributions need a density type array"
                 
                cmodel = ghmmwrapper.new_smodel()
                cmodel.M = len(B[0][0]) # Number of mixture componenents for emission distribution
                cmodel.prior = 0 # Unused
                
                # determining number of transition classes
                cos = 0
                if type(A[0][0]) == list:
                    cos = len(A)
                    cmodel.N = len(A[0])
                    # allocating class switching context
                    ghmmwrapper.ghmm_c_class_change_alloc(cmodel)
                    
                else: 
                    cos = 1
                    cmodel.N = len(A)
                    A = [A]
                    
                cmodel.cos = ghmmhelper.classNumber(A)  # number of transition classes in GHMM
                
                states = ghmmwrapper.arraysstate(cmodel.N)

                # Switching functions and transition classes are handled
                # elswhere

                #initialize states
                for i in range(cmodel.N):
                    state = ghmmwrapper.get_sstate_ptr(states,i)
                    state.pi = pi[i]
		    state.M = len(B[0][0])

                    # allocate arrays of emmission parameters
                    state.c = ghmmhelper.list2arrayd([1.0]) # Mixture weight
                    mu_list = B[i][0]
                    sigma_list = B[i][1]
                    a_list = B[i][2]
                    weight_list = B[i][3]
                    
                    state.mue = ghmmhelper.list2arrayd(mu_list) #mu = mue in GHMM C-lib.
                    state.u = ghmmhelper.list2arrayd(sigma_list)
                    state.c = ghmmhelper.list2arrayd(weight_list)
		    state.a = ghmmhelper.list2arrayd(a_list)

                    # setting densities types (all normal by default)                    
                    densit = ghmmwrapper.arraydensity(cmodel.M)
                    state.density = densit

                    mix_fix = [0] * state.M

                    for j in range(state.M):                   
                      ghmmwrapper.set_density(state,j,densities[i][j])
                      if densities[i][j] == ghmmwrapper.uniform:
                        mix_fix[j] = 1

                    # mixture fixing deactivated by default
                    state.mixture_fix = ghmmhelper.list2arrayint(mix_fix)
                    
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

class BackgroundDistribution:
    def __init__(self, emissionDomain, bgInput):
        
        if type(bgInput) == list:
            self.emissionDomain = emissionDomain
            distNum = len(bgInput)
        
            order = ghmmwrapper.int_array(distNum)
            b = ghmmwrapper.double_2d_array_nocols(distNum)
            for i in range(distNum):
                if len(emissionDomain) > 1:
                    o = math.log(len(bgInput[i]), len(emissionDomain)) - 1
                else:
                    o = len(bgInput[i]) - 1
                         
                assert (o % 1) == 0, "Invalid order of distribution " + str(i) + ": " + str(o)

                ghmmwrapper.set_arrayint(order, i, int(o))
                # dynamic allocation, rows have different lenghts
                b_i = ghmmhelper.list2arrayd(bgInput[i])
                ghmmwrapper.set_2d_arrayd_col(b,i,b_i)
    
            self.cbackground = ghmmwrapper.ghmm_d_background_alloc(distNum, len(emissionDomain), order, b)

        elif isinstance(bgInput, ghmmwrapper.background_distributions):
            self.cbackground = bgInput
            self.emissionDomain = emissionDomain
            
        else:
            raise TypeError, "Input type "+str(type(bgInput)) +" not recognized."    

    def __del__(self):
        log.debug( "__del__ BackgroundDistribution " + str(self.cbackground))
        ghmmwrapper.ghmm_d_background_free(self.cbackground)
        self.cbackground = None
    
    def __str__(self):
        outstr = "BackgroundDistribution instance:\n"
        outstr += "Number of distributions: " + str(self.cbackground.n)+"\n\n"
        outstr += str(self.emissionDomain) + "\n"
        d = ghmmhelper.matrixd2list(self.cbackground.b, self.cbackground.n, len(self.emissionDomain))
        outstr += "Distributions:\n"   
        for i in range(self.cbackground.n):
            outstr += "  Order: " + str(ghmmwrapper.get_arrayint(self.cbackground.order, i))+"\n"
            outstr += "  " + str(i+1) +": "+str(d[i])+"\n"
            
                     
        return outstr
        
    def tolist(self):
        dim = self.cbackground.m
        distNum = self.cbackground.n
        orders = ghmmhelper.arrayint2list(self.cbackground.order, distNum)
        B = []
        for i in xrange(distNum):
             order = orders[i]
             size = int(pow(m,(order+1)))
             b = [0.0]*size
             for j in xrange(size):
                  b[j] = ghmmwrapper.get_2d_arrayd(self.cbackground.b,i,j)
             B.append(b)
        return (distNum,orders,B)

#-------------------------------------------------------------------------------
#- HMM and derived  
class HMM:
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

        
        # Function pointers to the C level (assigned in derived classes)
        self.freeFunction = ""           # C function for deallocation
        self.samplingFunction = ""       # C function for sequence generation
        self.viterbiFunction = ""        # C function viterbi algorithm
        self.forwardFunction = ""        # C function forward algortihm (likelihood only)
        self.forwardAlphaFunction = ""   # C function forward algorithm (alpha matrix,scale)
        self.backwardBetaFunction = ""   # C function backkward algorithm (beta matrix)
        self.getStatePtr = ""            # C function to get a pointer to a state struct
        self.getModelPtr = ""            # C function to get a pointer to the model struct
        self.castModelPtr = ""           # C function to cast a *model to **model
        self.distanceFunction = ""       # C function to compute a probabilistic distance between models

    def __del__(self):
        """ Deallocation routine for the underlying C data structures. """
        log.debug( "__del__ HMM" + str(self.cmodel))
        self.freeFunction(self.cmodel)
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
        
        if not isinstance(emissionSequences,EmissionSequence) and not isinstance(emissionSequences,SequenceSet):
            raise TypeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)
        

        
        logPsum = sum(self.loglikelihoods(emissionSequences) )

        return logPsum

    def loglikelihoods(self, emissionSequences): 
        """ Compute a vector ( log( P[s| model]) )_{s} of log-likelihoods of the
            individual emission_sequences using the forward algorithm

            emission_sequences is of type SequenceSet

            Result: log( P[emissionSequences| model]) of type float
                    (numarray) vector of floats

        """

        if isinstance(emissionSequences,EmissionSequence):
            seqNumber = 1 
        elif isinstance(emissionSequences,SequenceSet):
            seqNumber = len(emissionSequences)        
        else:    
            raise TypeError, "EmissionSequence or SequenceSet required, got " + \
                  str(emissionSequences.__class__.__name__)        
              

        likelihood = ghmmwrapper.double_array(1)
        likelihoodList = []

        for i in range(seqNumber):
            seq = emissionSequences.getPtr(emissionSequences.cseq.seq,i)
            tmp = ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i)
            
            ret_val = self.forwardFunction(self.cmodel, seq, tmp, likelihood)
            if ret_val == -1:
                
                log.warning("forward returned -1: Sequence"+str(i)+"cannot be build.")
                # XXX Eventually this should trickle down to C-level
                # Returning -DBL_MIN instead of infinity is stupid, since the latter allows
                # to continue further computations with that inf, which causes
                # things to blow up later.
                # forwardFunction could do without a return value if -Inf is returned
                # What should be the semantics in case of computing the likelihood of
                # a set of sequences
                likelihoodList.append(-float('Inf'))
            else:
                likelihoodList.append(ghmmwrapper.get_arrayd(likelihood,0))

        ghmmwrapper.freearray(likelihood)
        likelihood = None
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
        pass
        

    def statePosterior(self, sequence, state, time):
        """ Return the log posterior probability for being at 'state' at time 'time' in 'sequence'.
            
            CAVEAT: statePosterior needs to calculate the complete forward and backward matrices. If 
            you are interested in multiple states it would be more efficient to use the posterior function
            directly and not multiple calls to statePosterior
        
        """
        # implemented in derived classes 
        pass


    def posterior(self, sequence):
        """ Posterior distribution matrix for 'sequence'.
        
        """
        # implemented in derived classes 
        pass



    def logprob(self, emissionSequence, stateSequence):
        """log P[ emissionSequence, stateSequence| m] 
        
            Defined in derived classes.
        """
        pass

    # The functions for model training are defined in the derived classes.
    def baumWelch(self, trainingSequences, nrSteps, loglikelihoodCutoff):
        pass

    def baumWelchSetup(self, trainingSequences, nrSteps):
        pass

    def baumWelchStep(self, nrSteps, loglikelihoodCutoff):
        pass
    
    def baumWelchDelete(self):
        pass
        
    # extern double ghmm_c_prob_distance(smodel *cm0, smodel *cm, int maxT, int symmetric, int verbose);
    def distance(self,model, seqLength):
        """ Returns the distance between 'self.cmodel' and 'model'.   """
        return self.distanceFunction(self.cmodel, model.cmodel, seqLength,0,0)


    
    def forward(self, emissionSequence):
        """

            Result: the (N x T)-matrix containing the forward-variables
                    and the scaling vector
        """
                      
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        unused = ghmmwrapper.double_array(1) # Dummy return value for forwardAlphaFunction    	
     
        t = len(emissionSequence)
        calpha = ghmmwrapper.double_2d_array (t, self.N)
        cscale = ghmmwrapper.double_array(t)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)

        error = self.forwardAlphaFunction(self.cmodel, seq,t, calpha, cscale, unused)
        if error == -1:
            log.error( "forward finished with -1: EmissionSequence cannot be build.")

        # translate alpha / scale to python lists 
        pyscale = ghmmhelper.arrayd2list(cscale, t)
        pyalpha = ghmmhelper.matrixd2list(calpha,t,self.N)
        
        # deallocation
        ghmmwrapper.freearray(unused)
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(calpha,t)
        unused = None
        csale = None
        calpha = None
        return (pyalpha,pyscale)
        

    def backward(self, emissionSequence, scalingVector):
        """

            Result: the (N x T)-matrix containing the backward-variables
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)
        
        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2arrayd(scalingVector)
        
        # alllocating beta matrix
        t = len(emissionSequence)
        cbeta = ghmmwrapper.double_2d_array(t, self.N)
        
        error = self.backwardBetaFunction(self.cmodel,seq,t,cbeta,cscale)
        if error == -1:
            log.error( "backward finished with -1: EmissionSequence cannot be build.")
            
        
        pybeta = ghmmhelper.matrixd2list(cbeta,t,self.N)

        # deallocation
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(cbeta,t)
        cscale = None
        cbeta = None
        return pybeta


    def viterbi(self, emissionSequences):
        """ Compute the Viterbi-path for each sequence in emissionSequences

            emission_sequences can either be a SequenceSet or an EmissionSequence

            Result: [q_0, ..., q_T] the viterbi-path of emission_sequences is an emmissionSequence
            object, [[q_0^0, ..., q_T^0], ..., [q_0^k, ..., q_T^k]} for a k-sequence
                    SequenceSet
        """
        
        if isinstance(emissionSequences,EmissionSequence):
            seqNumber = 1
        elif isinstance(emissionSequences,SequenceSet):
            seqNumber = len(emissionSequences)        
        else:    
            raise TypeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)


        log_p = ghmmwrapper.double_array(1)

        allLogs = []
        allPaths = []
        for i in range(seqNumber):
            seq = emissionSequences.getPtr(emissionSequences.cseq.seq,i)
            seq_len = ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i)

            if seq_len > 0:
                viterbiPath = self.viterbiFunction(self.cmodel,seq,seq_len,log_p)
            else:
                viterbiPath = None

            onePath = []
            
            # for model types without possible silent states
            # the length of the viterbi path is known
            if (self.cmodel.model_type &  4) == 0:    # check model_type for silent state flag        
                for j in range(seq_len):                
                    onePath.append(ghmmwrapper.get_arrayint(viterbiPath,j))
            
            # in the silent case we have to append as long as the path is
            # positive because the final path position is marked with a -1
            # in the following array entry.
            elif self.cmodel.model_type &  4 != 0:  # check model_type for silent state flag   
                
                for j in range( seq_len * self.N): # maximum length of a viterbi path for a silent model
                    d = ghmmwrapper.get_arrayint(viterbiPath,j)
                                   
                    if d >= 0:
                        onePath.append(d)
                    else:
                        break
                        
            allPaths.append(onePath)
            allLogs.append(ghmmwrapper.get_arrayd(log_p, 0))
            ghmmwrapper.freearray(viterbiPath) 
            viterbiPath = None

        
        ghmmwrapper.freearray(log_p)
        log_p = None
            
        if emissionSequences.cseq.seq_number > 1:
            return (allPaths, allLogs)
        else:
            return (allPaths[0], allLogs[0])


    def sample(self, seqNr ,T,seed = 0):
        """ Sample emission sequences 
                seqNr = number of sequences to be sampled
                T = length of each sequence
                seed = initialization value for rng, default 0 means 

        """
        seqPtr = self.samplingFunction(self.cmodel,seed,T,seqNr,self.N)
        seqPtr.state_labels = None
        seqPtr.state_labels_len = None

        return SequenceSet(self.emissionDomain,seqPtr)
        

    def sampleSingle(self, T, seed = 0):
        """ Sample a single emission sequence of length at most T.
            Returns a Sequence object.
        """
        seqPtr = self.samplingFunction(self.cmodel,seed,T,1,self.N)
        seqPtr.state_labels = None
        seqPtr.state_labels_len = None

        return EmissionSequence(self.emissionDomain,seqPtr)

    def state(self, stateLabel):
        """ Given a stateLabel return the integer index to the state 

            (state labels not yet implemented)
        """
        pass

    def getInitial(self, i):
        """ Accessor function for the initial probability \pi_i """
        state = self.getStatePtr(self.cmodel.s,i)
        return state.pi

    def setInitial(self, i, prob, fixProb=0):
        """ Accessor function for the initial probability \pi_i
            For 'fixProb' = 1 \pi will be rescaled to 1 with 'pi[i]' fixed to the
            arguement value of 'prob'.

         """

        state = self.getStatePtr(self.cmodel.s,i)
        old_pi = state.pi
        state.pi = prob

        # renormalizing pi, pi(i) is fixed on value 'prob'
        if fixProb == 1:
            coeff = (1.0 - old_pi) / prob
            for j in range(self.N):
                if i != j:
                    state = self.getStatePtr(self.cmodel.s,j)
                    p = state.pi
                    state.pi = p / coeff

    def getTransition(self, i, j):
        """ Accessor function for the transition a_ij """
        state = self.getStatePtr(self.cmodel.s,i)	

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."
        assert 0 <= j < self.N, "Index " + str(j) + " out of bounds."

        transition = None
        for i in range(state.out_states):
            stateId = ghmmwrapper.get_arrayint(state.out_id,i)
            if stateId == j:
                transition = ghmmwrapper.get_arrayd(state.out_a,i)
                break
        if transition:
            return transition
        else:
            return 0

    def setTransition(self, i, j, prob):
        """ Accessor function for the transition a_ij. """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."
        assert 0 <= j < self.N, "Index " + str(j) + " out of bounds."

        ghmmwrapper.ghmm_d_transition_set(self.cmodel, i, j, prob)


    def getEmission(self, i):
        """ Accessor function for the emission distribution parameters of state 'i'.

            For discrete models the distribution over the symbols is returned,
            for continuous models a matrix of the form
            [ [mu_1, sigma_1, weight_1] ... [mu_M, sigma_M, weight_M]  ] is returned.

        """
        if self.emissionDomain.CDataType == "int": # discrete emissions.
            state = self.getStatePtr(self.cmodel.s,i)
            if self.cmodel.model_type & 16:         #kHigherOrderEmissions
                order = ghmmwrapper.get_arrayint(self.cmodel.order, i)
                emissions = ghmmhelper.arrayd2list(state.b, self.M**(order+1))
            else:                
                emissions = ghmmhelper.arrayd2list(state.b, self.M)
            return emissions

        elif self.emissionDomain.CDataType == "double": # continuous emissions
            state = self.getStatePtr(self.cmodel.s,i)
            emParam = []
            for i in range(self.M):
                mixComp = []
                mixComp.append(ghmmwrapper.get_arrayd(state.mue,i) )
                mixComp.append(ghmmwrapper.get_arrayd(state.u,i) )
                mixComp.append(ghmmwrapper.get_arrayd(state.c,i) )
                emParam.append(mixComp)
            return emParam

    def setEmission(self, i, distributionParemters):
        """ Set the emission distribution parameters

            Defined in derived classes.
         """
        pass


    def toMatrices(self):
        "To be defined in derived classes."
        pass        


    def normalize(self):
        """ Normalize transition probs, emission probs (if applicable)

            Defined in derived classes.
        """
        pass


    def randomize(self, noiseLevel):
        """ """
        pass

    def write(self,fileName):
        """ Writes HMM to file 'fileName'.

        """
        self.fileWriteFunction(fileName,self.cmodel)

    def printtypes(self, model_type):
        strout = []
        first = 0
        if model_type &  2:         #kLeftRight
            strout.append("kLeftRight ")
        if model_type &  4:         #kSilentStates
            strout.append("kSilentStates ")
        if model_type &  8:         #kTiedEmissions
            strout.append("kTiedEmissions ")
        if model_type & 16:         #kHigherOrderEmissions
            strout.append("kHigherOrderEmissions ")
        if model_type & 32:         #kBackgroundDistributions
            strout.append("kBackgroundDistributions ")
        if model_type & 64:         #kClassLabels
            strout.append("kClassLabels ")
        if model_type == 0:         #kNotSpecified
            strout = "kNotSpecified"
        return join(strout,'')


def HMMwriteList(fileName,hmmList,fileType=GHMM_FILETYPE_SMO):
    if(fileType==GHMM_FILETYPE_SMO):
      if os.path.exists(fileName):
        log.warning( "HMMwriteList: File " + str(fileName) + " already exists. New models will be appended.")
      for model in hmmList:
        model.write(fileName)
    else:
       if os.path.exists(fileName):
         log.warning( "HMMwriteList: File " + str(fileName) + " already exists. Model will be overwritted.")
       models = ghmmwrapper.smodel_array(len(hmmList))
       for i,model in enumerate(hmmList):
         ghmmwrapper.set_smodel_ptr(models,model.cmodel,i)
       ghmmwrapper.ghmm_c_xml_write(fileName,models,len(hmmList))
       ghmmwrapper.free_smodel_array(models)
   
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
        
        # Assignment of the C function names to be used with this model type
        self.freeFunction = ghmmwrapper.call_ghmm_d_free
        self.samplingFunction = ghmmwrapper.ghmm_d_generate_sequences
        self.viterbiFunction = ghmmwrapper.ghmm_d_viterbi
        self.forwardFunction = ghmmwrapper.ghmm_d_forward_lean
        self.forwardAlphaFunction = ghmmwrapper.ghmm_d_forward      
        self.backwardBetaFunction = ghmmwrapper.ghmm_d_backward
        self.backwardTerminationFunction = ghmmwrapper.ghmm_d_backward_termination
        self.getStatePtr = ghmmwrapper.get_stateptr 
        self.fileWriteFunction = ghmmwrapper.call_model_print
        self.getModelPtr = ghmmwrapper.get_model_ptr
        self.castModelPtr = ghmmwrapper.cast_model_ptr
        self.distanceFunction = ghmmwrapper.ghmm_d_prob_distance
    
    def __del__(self):
        log.debug("__del__ DiscreteEmissionHMM" + str(self.cmodel))
        if self.cmodel.tied_to is not None:
            self.removeTiegroups()
        HMM.__del__(self)
        
    def __str__(self):
        hmm = self.cmodel
        strout = ["\nGHMM Model\n"]
        strout.append( "Name: " + str(self.cmodel.name))
        strout.append( "\nModelflags: "+ self.printtypes(self.cmodel.model_type))
        strout.append(  "\nNumber of states: "+ str(hmm.N))
        strout.append(  "\nSize of Alphabet: "+ str(hmm.M))
        if self.cmodel.model_type &  16: #kHigherOrderEmissions
            order = ghmmhelper.arrayint2list(self.cmodel.order, self.N)
        else:
            order = [0]*hmm.N

        for k in range(hmm.N):
            state = ghmmwrapper.get_stateptr(hmm.s, k)
            strout.append( "\n\nState number "+ str(k) +":")
            strout.append( "\nState order: " + str(order[k]))
            strout.append( "\nInitial probability: " + str(state.pi))
            #strout.append("\nsilent state: " + str(get_arrayint(self.cmodel.silent,k)))
            strout.append( "\nOutput probabilites: ")
            for outp in range(hmm.M**(order[k]+1)):
                strout.append(str(ghmmwrapper.get_arrayd(state.b,outp)))
                if outp % hmm.M == hmm.M-1:
                    strout.append( "\n")
                else:
                    strout.append( ", ")

            strout.append( "\nOutgoing transitions:")
            for i in range( state.out_states):
                strout.append( "\ntransition to state " + str( ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.out_a,i)))
            strout.append( "\nIngoing transitions:")
            for i in range(state.in_states):
                strout.append( "\ntransition from state " + str( ghmmwrapper.get_arrayint(state.in_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.in_a,i)))
            strout.append( "\nint fix:" + str(state.fix) + "\n")

        if self.cmodel.model_type &  4:
            strout.append("\nSilent states: \n")
            for k in range(hmm.N):
                strout.append( str(ghmmwrapper.get_arrayint(self.cmodel.silent,k)) + ", ")
        strout.append( "\n")
        return join(strout,'')


    def extendDurations(self, durationlist):
        """ extend states with durations larger one
            this done by explicit state copying in C """

        for i in range(len(durationlist)):
            if durationlist[i] > 1:
                error = ghmmwrapper.ghmm_d_duration_apply(self.cmodel, i, durationlist[i])
                self.N = self.cmodel.N
                if error:
                    log.error( "durations not applied")

    def setEmission(self, i, distributionParameters):
        """ Set the emission distribution parameters for a discrete model."""
        assert len(distributionParameters) == self.M
        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = self.getStatePtr(self.cmodel.s,i)

        # updating silent flag and/or model type if necessary 
        if self.cmodel.model_type & 4:
            if sum(distributionParameters) == 0.0:
                ghmmwrapper.set_arrayint(self.cmodel.silent, i, 1)
            else:                
                ghmmwrapper.set_arrayint(self.cmodel.silent, i, 0)
                #change model_type and free array if no silent state is left
                if 0 == sum(ghmmhelper.arrayint2list(self.cmodel.silent,self.N)):
                    self.cmodel.model_type -= 4
                    ghmmwrapper.freearray(self.cmodel.silent)
                    self.cmodel.silent = None
        #if the state becomes the first silent state allocate memory and set the silen flag
        elif sum(distributionParameters) == 0.0:
            self.cmodel.model_type += 4
            slist = [0]*self.N
            slist[i] = 1
            self.cmodel.silent = ghmmhelper.list2arrayint(slist)

        #set the emission probabilities
        for i in range(self.M):
            ghmmwrapper.set_arrayd(state.b,i,distributionParameters[i])


    def backwardTermination(self, emissionSequence, pybeta, scalingVector):
        """

            Result: the backward log probability of emissionSequence
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str (emissionSequence.__class__.__name__)

        seq = emissionSequence.getPtr (emissionSequence.cseq.seq, 0)
        
        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2arrayd (scalingVector)
        
        # alllocating beta matrix
        t = len (emissionSequence)
        cbeta = ghmmhelper.list2matrixd (pybeta)
        #print cbeta[0]

        # allocating double * for log probability
        log_p = ghmmwrapper.double_array (1)
        
        error = self.backwardTerminationFunction (self.cmodel, seq, t, cbeta[0], cscale, log_p)
        if error == -1:
            log.error("backward finished with -1: EmissionSequence cannot be build.")

        logp = ghmmwrapper.get_arrayd (log_p, 0)

        # deallocation
        ghmmwrapper.freearray(log_p)
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(cbeta[0],t)
        cscale = None
        cbeta = None
        return logp
    
    
    def logprob(self, emissionSequence, stateSequence):
        """ log P[ emissionSequence, stateSequence| m] """
        
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        state = self.getStatePtr( self.cmodel.s, stateSequence[0] )
        emissionProb = ghmmwrapper.get_arrayd(state.b, emissionSequence[0])
        if emissionProb == 0:
            silent = (self.cmodel.model_type & 4) and ghmmwrapper.get_arrayint(self.cmodel.silent, stateSequence[0])
            if silent == 1:
                emissionProb = 1
            else:
                raise SequenceCannotBeBuild, "first symbol " + str(emissionSequence[i+1]) + " not emitted by state " + str(stateSequence[0])
                        
        logP = math.log(state.pi * emissionProb )
        
        symbolIndex = 1

        try:
        
            for i in range(len(emissionSequence)-1):
                cur_state = self.getStatePtr( self.cmodel.s, stateSequence[i] )
                next_state = self.getStatePtr( self.cmodel.s, stateSequence[i+1] )
                for j in range(cur_state.out_states):
                    out_id = ghmmwrapper.get_arrayint(cur_state.out_id,j)
                    if out_id == stateSequence[i+1]:
                        emissionProb = ghmmwrapper.get_arrayd(next_state.b, emissionSequence[symbolIndex])
                        # print "b["+str(emissionSequence[symbolIndex])+"] in state " + str(stateSequence[i+1]) + " = ",emissionProb
                        symbolIndex += 1
                        if emissionProb == 0:
                            silent = (self.cmodel.model_type & 4) and ghmmwrapper.get_arrayint(self.cmodel.silent, stateSequence[i+1] )
                            if silent == 1:
                                emissionProb = 1
                                symbolIndex -= 1
                            else:
                                raise SequenceCannotBeBuild, "symbol " + str(emissionSequence[i+1]) + " not emitted by state "+ str(stateSequence[i+1])

                        logP += math.log( ghmmwrapper.get_arrayd(state.out_a,j) * emissionProb)
                        break
        except IndexError:
            pass
        return logP     
    
    def baumWelch(self, trainingSequences, nrSteps = None, loglikelihoodCutoff = None):
        """ Reestimates the model with the sequence in 'trainingSequences'.
           
            Note that training for models including silent states is not yet supported.

            nrSteps is the maximal number of BW-steps
            loglikelihoodCutoff is the least relative improvement in likelihood with respect to the last iteration 
            required to continue.
           
        """
        if not isinstance(trainingSequences,EmissionSequence) and not isinstance(trainingSequences,SequenceSet):
            raise TypeError, "EmissionSequence or SequenceSet required, got " + str(trainingSequences.__class__.__name__)

        if self.cmodel.model_type & 4:     #kSilentStates
            log.critical( "Sorry, training of models containing silent states not yet supported.")
        else:
            if nrSteps == None:
                ghmmwrapper.ghmm_d_baum_welch(self.cmodel, trainingSequences.cseq)
            else:
                assert loglikelihoodCutoff != None
                ghmmwrapper.ghmm_d_baum_welch_nstep(self.cmodel, trainingSequences.cseq,
                                                        nrSteps, loglikelihoodCutoff)


    def applyBackground(self, backgroundWeight):
        """Apply the background distribution to the emission probabilities of states which
           have been assigned one (usually in the editor and coded in the XML).
           applyBackground computes a convex combination of the emission probability and
           the background, where the backgroundWeight parameter (within [0,1]) controls
           the background's contribution for each state.
        """
        assert len(backgroundWeight) == self.N, "Argument 'backgroundWeight' does not match number of states."
        
        cweights = ghmmhelper.list2arrayd(backgroundWeight)
        result = ghmmwrapper.ghmm_d_background_apply(self.cmodel, cweights)
        
        ghmmwrapper.freearray(cweights)
        if result is not 0:
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
            ghmmwrapper.freearray(self.cmodel.background_id)
        self.cmodel.bp = backgroundObject.cbackground
        self.background = backgroundObject
        self.cmodel.background_id = ghmmhelper.list2arrayint(stateBackground)

        # updating model type
        self.cmodel.model_type += 32 #kBackgroundDistributions
    
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
            ghmmwrapper.set_arrayint(self.cmodel.background_id,i,stateBackground[i])
        
    def assignStateBackground(self, state, backgroundID):

        assert self.cmodel.background_id is not None, "Error: No backgrounds defined in model."   
        if self.labelDomain.isAdmissable(backgroundID):
            ghmmwrapper.set_arrayint(self.cmodel.background_id, state, backgroundID)
        else:
            log.error( str(backgroundID) + " is not contained in labelDomain."  )
    
    
    def getBackgroundAssignments(self):
        return ghmmhelper.arrayint2list(self.cmodel.background_id, self.N)
        
    
    def updateTieGroups(self):
        
        assert self.cmodel.tied_to is not None, "cmodel.tied_to is undefined."
        ghmmwrapper.ghmm_d_update_tied_groups(self.cmodel)

    
    def setTieGroups(self, tieList):
        
        assert len(tieList) == self.N, "Number of entries in tieList is different from number of states."
        
        if self.cmodel.tied_to is None:
            log.debug( "allocating tied_to")
            self.cmodel.tied_to = ghmmhelper.list2arrayint(tieList)
            self.cmodel.model_type += 8
        else:
            log.debug( "tied_to already initialized")
            for i in range(self.N):
                ghmmwrapper.set_arrayint(self.cmodel.tied_to,i,tieList[i])

    def removeTiegroups(self):
        ghmmwrapper.freearray(self.cmodel.tied_to)
        self.cmodel.tied_to = None
        self.cmodel.model_type -= 8
    
    def getTieGroups(self):
        assert self.cmodel.tied_to is not None, "cmodel.tied_to is undefined."
        
        return ghmmhelper.arrayint2list(self.cmodel.tied_to, self.N)
    
    def getSilentFlag(self,state):
        if self.cmodel.model_type & 4:
            return ghmmwrapper.get_arrayint(self.cmodel.silent,state)
        else:
            return 0
    

    def normalize(self):
        """ Normalize transition probs, emission probs (if applicable) """
        
        log.debug( "Normalizing now.")

        for i in range(self.N):
            # normalizing transitions
            state = self.getStatePtr(self.cmodel.s,i)
            pSum = 0.0
            stateIds = []
            for j in range(state.out_states):
                stateIds.append(ghmmwrapper.get_arrayint(state.out_id,j))
                pSum += ghmmwrapper.get_arrayd(state.out_a,j)
            for j in range(state.out_states):
                normP = ghmmwrapper.get_arrayd(state.out_a,j) / pSum
                ghmmwrapper.set_arrayd(state.out_a,j,normP) # updating out probabilities
                
                inState = self.getStatePtr(self.cmodel.s,stateIds[j])
                for k in range(inState.in_states):
                    inId = ghmmwrapper.get_arrayint(inState.in_id,k)
                    if inId == i:
                        ghmmwrapper.set_arrayd(inState.in_a,k,normP) # updating in probabilities

            # normalizing emissions    
            pSum = 0.0
            for j in range(self.M):
                pSum += ghmmwrapper.get_arrayd(state.b,j)

            for j in range(self.M):
                if pSum >0:  # check for silent state
                    normP = ghmmwrapper.get_arrayd(state.b,j) / pSum
                    ghmmwrapper.set_arrayd(state.b,j,normP)

    def toMatrices(self):
        "Return the parameters in matrix form."
        A = []
        B = []
        pi = []
        if self.cmodel.model_type & 16:         #kHigherOrderEmissions
            order = ghmmhelper.arrayint2list(self.cmodel.order, self.N)
        else:
            order = [0]*self.N

        for i in range(self.cmodel.N):
            A.append([0.0] * self.N)
            state = self.getStatePtr(self.cmodel.s,i)
            pi.append(state.pi)
            B.append(ghmmhelper.arrayd2list(state.b,self.M ** (order[i]+1)))
            for j in range(state.out_states):
                state_index = ghmmwrapper.get_arrayint(state.out_id,j)
                A[i][state_index] = ghmmwrapper.get_arrayd(state.out_a,j)

        return [A,B,pi]

    def isSilent(self,state):
        """ Returns True if 'state' is silent, False otherwise
        
        """
        assert 0 <= state <= self.N-1, "Invalid state index"
        
        # check whether model contains silent states at all
        if not (self.cmodel.model_type & 4):
            return False
        
        # check silent flag for state 'state'
        if ghmmwrapper.get_arrayint(self.cmodel.silent,state):
            return True
        else:
            return False    


    def pathPosterior(self, sequence, path):
        """ Returns the log posterior probability for 'path' having generated 'sequence'.
            
            CAVEAT: statePosterior needs to calculate the complete forward and backward matrices. If 
            you are interested in multiple paths it would be more efficient to use the 'posterior' function
            directly and not multiple calls to pathPosterior
        """
        # XXX for silent states things arr more complicated -> to be done
        if self.cmodel.model_type & 4:
            raise RuntimeError, "Models with silent states not yet supported."
            
                    
        # checking function arguments
        assert isinstance(sequence, EmissionSequence), "Input to posterior must be EmissionSequence object"
        

        # checking path validity (XXX too inefficient ?)
        for p in path:
            assert 0 <= p <= self.N-1, "Invalid state index "+str(p)+". Model and path are incompatible"

        # calculate complete posterior matrix
        post = self.posterior(sequence)
        path_posterior = []
        
        if not (self.cmodel.model_type & 4):   # check for kSilentStates
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
            you are interested in multiple states it would be more efficient to use the posterior function
            directly and not multiple calls to statePosterior
        
        """
        # XXX for silent states things arr more complicated -> to be done
        if self.cmodel.model_type & 4:
            raise RuntimeError, "Models with silent states not yet supported."
            
                    
        # checking function arguments
        assert isinstance(sequence, EmissionSequence), "Input to posterior must be EmissionSequence object"
        assert 0 <= time <= len(sequence), "Invalid sequence index: "+str(time)+" (sequence has length "+str(len(sequence))+" )."
        assert 0 <= state <= self.N-1, "Invalid state index: " +str(state)+ " (models has "+str(self.N)+" states )."

        post = self.posterior(sequence)
        return post[time][state]



    def posterior(self, sequence):
        """ Posterior distribution matrix for 'sequence'.
        
        """

        # XXX for silent states things arr more complicated -> to be done    	
        if self.cmodel.model_type & 4:
            raise RuntimeError, "Models with silent states not yet supported."

        
        assert isinstance(sequence, EmissionSequence), "Input to posterior must be EmissionSequence object"
        
        (alpha,scale)  = self.forward(sequence)
        beta = self.backward(sequence,scale)
       
#        print "alpha: "
#        for a in alpha:
#            print a
#        print "\n"
#        
#        print "beta: "
#        for b in beta:
#            print b
#        print "\n"    
        
        
        post_mat = []
        for i in range(len(sequence)):
            post = []            
            for state in range(self.N):
                #s = sum(alpha[i])
                #print alpha[i]
                post.append( alpha[i][state] * beta[i][state] ) 
       
            post_mat.append(post)   

       

        return post_mat


    def toXML(self, filename, backgroundobj = None):
        [A,B,pi] = self.toMatrices()
        nalpha = self.cmodel.M
        nstates = self.cmodel.N

        hmm_dom = xmlutil.HMM()

        hmm_dom.hmmClass.addCode(0, "C1", xmlutil.ValidatingString("Switching class"))

        alphabet = self.emissionDomain
        # adapted to the pair HMM xmlutil with multiple alphabets
        hmm_dom.hmmAlphabets[0] = xmlutil.DiscreteHMMAlphabet()
        for c in alphabet.listOfCharacters:
            hmm_dom.hmmAlphabets[0].addCode(alphabet.index[c], c)

        if backgroundobj != None:
             (n,orders,b) = backgroundobj.tolist()
             background_id = self.getBackgroundAssignments()
             for i in xrange(n):
                  hmm_dom.backgroundDistributions.addDistribution(str(i), orders[i], b[i])

        try:
            tiedlist = self.getTieGroups()
        except AssertionError:
            log.warning( "Ignore tied groups")
            log.warning( "self.cmodel.tied_to not defined")

        for i in xrange(self.cmodel.N):
            cstate = self.getStatePtr(self.cmodel.s,i)
            if nalpha > 1:
                order = int(math.log(len(B[i]), nalpha)) - 1
            else:
                order = len(B[i]) - 1

            have_background = False # XXX
            tiedto = None           # XXX
           
            state_dom = xmlutil.HMMState(-1,hmm_dom)
            state_dom.fromDiscreteState( i, pi[i], B[i], 0, order, tiedto, have_background)
            hmm_dom.state[state_dom.index] = state_dom
            hmm_dom.id2index[i] = state_dom.index
            
        nr_classes = 1
        hmm_dom.G.edgeWeights[0] = EdgeWeight(hmm_dom.G)

        #
        # Add transition probabilities
        #
        for i in xrange(self.cmodel.N):
            v = hmm_dom.id2index[i]
            for j in  xrange(self.cmodel.N):
                if A[i][j] > 0.0:
                    w = hmm_dom.id2index[j]
                    hmm_dom.G.AddEdge(v,w,A[i][j])

        hmm_dom.WriteXML(filename)
        del hmm_dom

######################################################
        
class StateLabelHMM(DiscreteEmissionHMM):
    """ Labelled HMMs with discrete emissions. 
        Same feature list as in DiscreteEmission models.    
    
    """
    def __init__(self, emissionDomain, distribution, labelDomain, cmodel):
        DiscreteEmissionHMM.__init__(self, emissionDomain, distribution, cmodel)
        assert isinstance(labelDomain, LabelDomain), "Invalid labelDomain"
        self.labelDomain = labelDomain
        
        # Assignment of the C function names to be used with this model type
        self.forwardFunction = ghmmwrapper.ghmm_d_logp
        self.forwardAlphaLabelFunction = ghmmwrapper.ghmm_dl_forward
        self.backwardBetaLabelFunction = ghmmwrapper.ghmm_dl_backward
        self.kbestFunction = ghmmwrapper.ghmm_dl_kbest        
        self.gradientDescentFunction = ghmmwrapper.ghmm_dl_gradient_descent
        self.cmodel.label = ghmmwrapper.int_array(self.N)

    def __str__(self):
        hmm = self.cmodel
        strout = ["\nGHMM Model\n"]
        strout.append("Name: " + str(self.cmodel.name))
        strout.append("\nModelflags: "+ self.printtypes(self.cmodel.model_type))
        strout.append("\nNumber of states: "+ str(hmm.N))
        strout.append("\nSize of Alphabet: "+ str(hmm.M))
        order = ghmmhelper.arrayint2list(hmm.order, self.N)
        label = ghmmhelper.arrayint2list(hmm.label, self.N)
        for k in range(hmm.N):
            state = ghmmwrapper.get_stateptr(hmm.s, k)
            strout.append("\n\nState number "+ str(k) +":")

            strout.append("\nState label: "+str(self.labelDomain.external(label[k])))

            strout.append("\nState order: " + str(order[k]))
            strout.append("\nInitial probability: " + str(state.pi))
            strout.append("\nOutput probabilites:\n")
            for outp in range(hmm.M**(order[k]+1)):
                strout+=str(ghmmwrapper.get_arrayd(state.b,outp))
                if outp % hmm.M == hmm.M-1:
                    strout.append("\n")
                else:
                    strout.append(", ")

            strout.append("Outgoing transitions:")
            for i in range( state.out_states):
                strout.append("\ntransition to state " + str( ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.out_a,i)))
            strout.append( "\nIngoing transitions:")
            for i in range(state.in_states):
                strout.append( "\ntransition from state " + str( ghmmwrapper.get_arrayint(state.in_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.in_a,i)))
                strout.append("\nint fix:" + str(state.fix) + "\n")

        if hmm.model_type & 4:
            strout.append("\nSilent states: \n")
            for k in range(hmm.N):
                strout.append(str(ghmmwrapper.get_arrayint(hmm.silent,k)) + ", ")
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
            
        ghmmwrapper.freearray(self.cmodel.label)
        self.cmodel.label = ghmmhelper.list2arrayint([self.labelDomain.internal(l) for l in labelList])

    def getLabels(self):
        labels = ghmmhelper.arrayint2list(self.cmodel.label, self.N)
        return [self.labelDomain.external(l) for l in labels]
    
    def getLabel(self,stateIndex):
        """ Returns label of the state 'stateIndex'.
        
        """
        return ghmmwrapper.get_arrayint(self.cmodel.label, stateIndex)
         
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

    def sampleSingle(self, seqLength,seed = 0):
        seqPtr = ghmmwrapper.ghmm_dl_generate_sequences(self.cmodel,seed,seqLength,1,seqLength)
        return EmissionSequence(self.emissionDomain, seqPtr, labelDomain = self.labelDomain )

    def sample(self, seqNr,seqLength, seed = 0):
        seqPtr = ghmmwrapper.ghmm_dl_generate_sequences(self.cmodel,seed,seqLength,seqNr,seqLength)
        return SequenceSet(self.emissionDomain,seqPtr, labelDomain = self.labelDomain)


    def viterbiLabels(self, emissionSequences):
        """ Returns on approximation of the most likely labeling of the input sequence(s) as
            given by the viterbi path.
        """

        if isinstance(emissionSequences,EmissionSequence):
            seqNumber = 1
        elif isinstance(emissionSequences,SequenceSet):
            seqNumber = len(emissionSequences)
        else:
            raise TypeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)

        assert emissionSequences.emissionDomain == self.emissionDomain, "Sequence and model emissionDomains are incompatible."
        
        (vPath, log_p) = self.viterbi(emissionSequences)

        labels = []

        if seqNumber == 1:
            labels = map(lambda i: self.labelDomain.external(self.getLabel(i)), vPath)
        else:
            for j in range(seqNumber):
                labels.append( map(lambda i: self.labelDomain.external(self.getLabel(i)), vPath[j]) )

        return (labels, log_p)


    def kbest(self, emissionSequences, k = 1):
        """ Compute the k probable labeling for each sequence in emissionSequences

            emissionSequences can either be a SequenceSet or an EmissionSequence

            Result: [l_0, ..., l_T] the labeling of emissionSequences is an emmissionSequence
            object, [[l_0^0, ..., l_T^0], ..., [l_0^k, ..., l_T^k]} for a k-sequence
                    SequenceSet
        """
        if self.cmodel.model_type & 4:     #kSilentStates
            log.critical( "Sorry, k-best decoding on models containing silent states not yet supported.")
        else:
            if isinstance(emissionSequences,EmissionSequence):
                seqNumber = 1
            elif isinstance(emissionSequences,SequenceSet):
                seqNumber = len(emissionSequences)
            else:
                raise TypeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)

            log_p = ghmmwrapper.double_array(1)

            allLogs = []
            allLabels = []

            for i in range(seqNumber):
                seq = emissionSequences.getPtr(emissionSequences.cseq.seq, i)
                seq_len = ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i)

                labeling =  self.kbestFunction(self.cmodel, seq, seq_len, k, log_p)
                oneLabel = []

                for j in range(seq_len):
                    oneLabel.append(ghmmwrapper.get_arrayint(labeling,j))

                allLabels.append(oneLabel)
                allLogs.append(ghmmwrapper.get_arrayd(log_p, 0))
                ghmmwrapper.freearray(labeling)

                            
            ghmmwrapper.freearray(log_p)
            labeling = None
            log_p = None

            if emissionSequences.cseq.seq_number > 1:
                return (map(self.externalLabel, allLabels), allLogs)
            else:
                return (self.externalLabel(allLabels[0]), allLogs[0])


    def gradientSearch(self, emissionsequences, eta=.1, steps=20):
        """ trains a model with given sequencesgradescentFunction using gradient descent

            emission_sequences can either be a SequenceSet or an EmissionSequence

        """
        
        # check for labels 
        assert self.cmodel.model_type & 64, "Error: Model is not labelled. "
        
        
        if isinstance(emissionsequences, EmissionSequence):
            seqNumber = 1
        elif isinstance(emissionsequences, SequenceSet):
            seqNumber = len(emissionsequences)
        else:
            raise TypeError, "LabeledEmissionSequence or LabeledSequenceSet required, got "\
                  + str(emissionsequences.__class__.__name__)


        cmodelPTR = ghmmwrapper.modelarray_alloc(1)
        ghmmwrapper.modelarray_setptr(cmodelPTR, self.cmodel, 0)
        error = self.gradientDescentFunction(cmodelPTR, emissionsequences.cseq, eta, steps)

        self.cmodel = ghmmwrapper.modelarray_getptr(cmodelPTR, 0)
        ghmmwrapper.freearray(cmodelPTR)
        
        if error == -1:
            log.error("Gradient descent finished not successfully.")

        return error
        
    def modelNormalize(self):
       i_error = ghmmwrapper.ghmm_d_normalize(self.cmodel)
       if i_error == -1:
            log.error("normalize finished with -1")
       

        
    def labelSeqLikelihoods(self, emissionSequences):
        """ Compute a vector ( log( P[s,l| model]) )_{s} of log-likelihoods of the
            individual emission_sequences using the forward algorithm

            emission_sequences is of type SequenceSet

            Result: log( P[emissionSequences,labels| model]) of type float
                    (numarray) vector of floats

        """
        if isinstance(emissionSequences,EmissionSequence):
            seqNumber = 1
        elif isinstance(emissionSequences,SequenceSet):
            seqNumber = len(emissionSequences)
        else:
            raise TypeError, "EmissionSequence or SequenceSet required, got " + \
                  str(emissionSequences.__class__.__name__)

        assert emissionSequences.cseq.state_labels is not None, "Sequence needs to be labeled."

        likelihood = ghmmwrapper.double_array(1)
        likelihoodList = []


        for i in range(seqNumber):
            seq = emissionSequences.getPtr(emissionSequences.cseq.seq,i)
            labels = ghmmwrapper.get_col_pointer_int(emissionSequences.cseq.state_labels,i)
            tmp = ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i)
            ret_val = ghmmwrapper.ghmm_dl_logp(self.cmodel, seq, labels, tmp, likelihood)

            if ret_val == -1:

                log.warning("forward returned -1: Sequence"+ str(i) +"cannot be build.")
                likelihoodList.append(-float('Inf'))
            else:
                likelihoodList.append(ghmmwrapper.get_arrayd(likelihood,0))

        ghmmwrapper.freearray(likelihood)
        likelihood = None
        return likelihoodList

    def forwardLabels(self, emissionSequence, labelSequence):
        """

            Result: the (N x T)-matrix containing the forward-variables
                    and the scaling vector
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        logP = ghmmwrapper.double_array(1)
        n_states = self.cmodel.N

        t = ghmmwrapper.get_arrayint(emissionSequence.cseq.seq_len,0)
        if t != len(labelSequence):
            raise TypeError, "ERROR: Observation and Labellist must have same length"

        calpha = ghmmwrapper.double_2d_array(t, n_states)
        cscale = ghmmwrapper.double_array(t)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)
        label = ghmmwrapper.int_array(t)

        for i in range(len(labelSequence)):
            ghmmwrapper.set_arrayint(label, i, self.internalLabel(labelSequence[i]))

        error = self.forwardAlphaLabelFunction(self.cmodel, seq, label, t, calpha, cscale, logP)
        if error == -1:
            log.error( "Forward finished with -1: Sequence " + str(i) + " cannot be build.")

        # translate alpha / scale to python lists
        pyscale = ghmmhelper.arrayd2list(cscale, t)
        pyalpha = ghmmhelper.matrixd2list(calpha,t,n_states)
        logpval = ghmmwrapper.get_arrayd(logP, 0)
        
        ghmmwrapper.freearray(label)
        ghmmwrapper.freearray(logP)
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(calpha,t)
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

        t = ghmmwrapper.get_arrayint(emissionSequence.cseq.seq_len,0)
        if t != len(labelSequence):
            raise TypeError, "ERROR: Observation and Labellist must have same length"

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)
        label = ghmmwrapper.int_array(t)

        for i in range(len(labelSequence)):
            ghmmwrapper.set_arrayint(label, i, self.internalLabel(labelSequence[i]))

        logP = ghmmwrapper.double_array(1)

        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2arrayd(scalingVector)

        # alllocating beta matrix
        cbeta = ghmmwrapper.double_2d_array(t, self.cmodel.N)

        error = self.backwardBetaLabelFunction(self.cmodel,seq,label,t,cbeta,cscale,logP)
        if error == -1:
            log.error( "backward finished with -1: EmissionSequence cannot be build.")
            

        pybeta = ghmmhelper.matrixd2list(cbeta,t,self.cmodel.N)
        logpval = ghmmwrapper.get_arrayd(logP, 0)
        
        # deallocation
        ghmmwrapper.freearray(logP)
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.freearray(label)
        ghmmwrapper.free_2darrayd(cbeta,t)
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

        if self.cmodel.model_type & 4:     #kSilentStates
            log.critical("Sorry, training of models containing silent states not yet supported.")
        else:
            if nrSteps == None:
                ghmmwrapper.ghmm_dl_baum_welch(self.cmodel, trainingSequences.cseq)
            else:
                assert loglikelihoodCutoff != None
                ghmmwrapper.ghmm_dl_baum_welch_nstep(self.cmodel, trainingSequences.cseq,
                                                        nrSteps, loglikelihoodCutoff)


class GaussianEmissionHMM(HMM):
    """ HMMs with Gaussian distribution as emissions.
        
    """

    def __init__(self, emissionDomain, distribution, cmodel):
        HMM.__init__(self, emissionDomain, distribution, cmodel)

        # Assignment of the C function names to be used with this model type
        self.freeFunction = ghmmwrapper.call_smodel_free
        self.samplingFunction = ghmmwrapper.ghmm_c_generate_sequences
        self.viterbiFunction = ghmmwrapper.ghmm_c_viterbi
        self.forwardFunction = ghmmwrapper.ghmm_c_logp
        self.forwardAlphaFunction = ghmmwrapper.ghmm_c_forward
        self.backwardBetaFunction = ghmmwrapper.ghmm_c_backward
        self.getStatePtr = ghmmwrapper.get_sstate_ptr
        self.fileWriteFunction = ghmmwrapper.call_smodel_print
        self.getModelPtr = ghmmwrapper.get_smodel_ptr
        self.castModelPtr = ghmmwrapper.cast_smodel_ptr
        self.distanceFunction = ghmmwrapper.ghmm_c_prob_distance

        # Baum Welch context, call baumWelchSetup to initalize
        self.BWcontext = ""

    def getTransition(self, i, j):
        """ Returns the probability of the transition from state i to state j.
             Raises IndexError if the transition is not allowed
        """
        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."
        assert 0 <= j < self.N, "Index " + str(j) + " out of bounds."
 
     	transition = ghmmwrapper.smodel_get_transition(self.cmodel, i, j, 0)
        if transition < 0.0: # Tried to access non-existing edge:
            transition = 0.0
        return transition

    def setTransition(self, i, j, prob):
        """ Accessor function for the transition a_ij """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."
        assert 0 <= j < self.N, "Index " + str(j) + " out of bounds."

        ghmmwrapper.smodel_set_transition(self.cmodel, i, j, 0, float(prob))

    def getEmission(self, i):
        """ Return (mu, sigma^2)  """

        # ensure proper index
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = ghmmwrapper.get_sstate(self.cmodel, i)
        mu = ghmmwrapper.get_arrayd(state.mue, 0)
        sigma = ghmmwrapper.get_arrayd(state.u,0)
        return (mu, sigma)

    def setEmission(self, i, (mu, sigma)):
        """ Set the emission distributionParameters for state i """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = ghmmwrapper.get_sstate(self.cmodel, i)
        ghmmwrapper.set_arrayd(state.mue, 0, float(mu))  # GHMM C is german: mue instead of mu
        sigma = ghmmwrapper.set_arrayd(state.u, 0, float(sigma))


    def getMixtureFix(self,state):
        s = self.getStatePtr(self.cmodel.s,state)
        return ghmmhelper.arrayint2list(s.mixture_fix,self.M)
        
        
    def setMixtureFix(self, state ,flags):    
        s = self.getStatePtr(self.cmodel.s,state)        
        s.mixture_fix = ghmmhelper.list2arrayint(flags)
    
    def getStateFix(self,state):
        s = self.getStatePtr(self.cmodel.s,state)
        return s.fix
        
    def setStateFix(self, state ,flag):    
        s = self.getStatePtr(self.cmodel.s,state)
        s.fix = flag
    
    def __str__(self):
        hmm = self.cmodel
        # TTTTTTTTTTTT XXX
        strout = ["\nHMM Overview:"]
        strout.append("\nNumber of states: " + str(hmm.N))
        strout.append("\nNumber of mixture components: " + str(hmm.M))

        for k in range(hmm.N):
            state = ghmmwrapper.get_sstate(hmm, k)
            strout.append("\n\nState number "+ str(k) + ":")
            strout.append("\nInitial probability: " + str(state.pi) + "\n")

            weight = ""
            mue = ""
            u =  ""

            weight += str(ghmmwrapper.get_arrayd(state.c,0))
            mue += str(ghmmwrapper.get_arrayd(state.mue,0))
            u += str(ghmmwrapper.get_arrayd(state.u,0))

            strout.append("  mean: " + str(mue) + "\n")
            strout.append("  variance: " + str(u) + "\n")
            strout.append("  fix: " + str(state.fix) + "\n")
            
            for c in range(self.cmodel.cos):
                strout.append("\n  Class : " + str(c)                )
                strout.append("\n    Outgoing transitions:")
                for i in range( state.out_states):
                    strout.append("\n      transition to state " + str(ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability = " + str(ghmmwrapper.get_2d_arrayd(state.out_a,c,i)))
                strout.append("\n    Ingoing transitions:")
                for i in range(state.in_states):
                    strout.append("\n      transition from state " + str(ghmmwrapper.get_arrayint(state.in_id,i) ) +" with probability = "+ str(ghmmwrapper.get_2d_arrayd(state.in_a,c,i)))


        return join(strout,'')

    # different function signatures require overloading of parent class methods    
    def sample(self, seqNr ,T,seed = -1):
        """ Sample emission sequences 
        
        """
        seqPtr = self.samplingFunction(self.cmodel,seed,T,seqNr,0,-1) 

        return SequenceSet(self.emissionDomain,seqPtr)


    def sampleSingle(self, T,seed= -1):
        """ Sample a single emission sequence of length at most T.
            Returns a Sequence object.
        """
        seqPtr = self.samplingFunction(self.cmodel,seed,T,1,0,-1) 

        return EmissionSequence(self.emissionDomain,seqPtr)

    def forward(self, emissionSequence):
        """

            Result: the (N x T)-matrix containing the forward-variables
                    and the scaling vector
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        logP = ghmmwrapper.double_array(1)
        i = self.cmodel.N


        t = ghmmwrapper.get_arrayint(emissionSequence.cseq.seq_len,0)
        calpha = ghmmwrapper.double_2d_array (t, i)
        cscale = ghmmwrapper.double_array(t)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)

        error = self.forwardAlphaFunction(self.cmodel, seq,t, None, calpha, cscale, logP)
        if error == -1:
            log.error( "Forward finished with -1: Sequence " + str(seq_nr) + " cannot be build.")
        
        # translate alpha / scale to python lists
        pyscale = ghmmhelper.arrayd2list(cscale, t)
        pyalpha = ghmmhelper.matrixd2list(calpha,t,i)

        ghmmwrapper.freearray(logP)
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(calpha,t)
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

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)

        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2arrayd(scalingVector)

        # alllocating beta matrix
        t = ghmmwrapper.get_arrayint(emissionSequence.cseq.seq_len,0)
        cbeta = ghmmwrapper.double_2d_array(t, self.cmodel.N)

        error = self.backwardBetaFunction(self.cmodel,seq,t,None,cbeta,cscale)
        if error == -1:
            log.error( "backward finished with -1: EmissionSequence cannot be build.")
            

        pybeta = ghmmhelper.matrixd2list(cbeta,t,self.cmodel.N)

        # deallocation
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(cbeta,t)
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

        if isinstance(emissionSequences,EmissionSequence):
            seqNumber = 1 
        elif isinstance(emissionSequences,SequenceSet):
            seqNumber = len(emissionSequences)        
        else:    
            raise TypeError, "EmissionSequence or SequenceSet required, got " + \
                  str(emissionSequences.__class__.__name__)        

        if self.cmodel.cos > 1:
            log.debug( "self.cmodel.cos = " + str( self.cmodel.cos) )
            assert self.cmodel.class_change is not None, "Error: class_change not initialized."

        likelihood = ghmmwrapper.double_array(1)
        likelihoodList = []

        for i in range(seqNumber):
            seq = emissionSequences.getPtr(emissionSequences.cseq.seq,i)
            tmp = ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i)
            
            if self.cmodel.cos > 1:
                self.cmodel.class_change.k = i
            
            ret_val = self.forwardFunction(self.cmodel, seq, tmp, likelihood)
            if ret_val == -1:
                
                log.warning( "forward returned -1: Sequence"+str(i)+"cannot be build.")
                # XXX Eventually this should trickle down to C-level
                # Returning -DBL_MIN instead of infinity is stupid, since the latter allows
                # to continue further computations with that inf, which causes
                # things to blow up later.
                # forwardFunction could do without a return value if -Inf is returned
                # What should be the semantics in case of computing the likelihood of
                # a set of sequences
                likelihoodList.append(-float('Inf'))
            else:
                likelihoodList.append(ghmmwrapper.get_arrayd(likelihood,0))

        ghmmwrapper.freearray(likelihood)
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
        
        if isinstance(emissionSequences,EmissionSequence):
            seqNumber = 1
        elif isinstance(emissionSequences,SequenceSet):
            seqNumber = len(emissionSequences)        
        else:    
            raise TypeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)

        if self.cmodel.cos > 1:
            log.debug( "self.cmodel.cos = "+ str( self.cmodel.cos))
            assert self.cmodel.class_change is not None, "Error: class_change not initialized."


        log_p = ghmmwrapper.double_array(1)

        allLogs = []
        allPaths = []
        for i in range(seqNumber):
            if self.cmodel.cos > 1:
                # if emissionSequence is a sequenceSet with multiple sequences, 
                # use sequence index as class_change.k
                if emissionSequences.cseq.seq_number > 1:
                    self.cmodel.class_change.k = i

            seq = emissionSequences.getPtr(emissionSequences.cseq.seq,i)
            seq_len = ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i)

            if seq_len != 0:
                viterbiPath = self.viterbiFunction(self.cmodel,seq,seq_len,log_p)
            else:
                viterbiPath = None

            onePath = []
            for j in range(seq_len):                
                onePath.append(ghmmwrapper.get_arrayint(viterbiPath,j))
            

            allPaths.append(onePath)
            allLogs.append(ghmmwrapper.get_arrayd(log_p, 0))
        
        ghmmwrapper.freearray(log_p)
        ghmmwrapper.freearray(viterbiPath)  
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

        for i in range(self.N):
            # normalizing transitions
            state = self.getStatePtr(self.cmodel.s,i)
            pSum = 0.0
            stateIds = []
            for j in range(state.out_states):
                stateIds.append(ghmmwrapper.get_arrayint(state.out_id,j))
                pSum += ghmmwrapper.get_arrayd(state.out_a,j)
            for j in range(state.out_states):
                normP = ghmmwrapper.get_arrayd(state.out_a,j) / pSum
                ghmmwrapper.set_arrayd(state.out_a,j,normP) # updating out probabilities

                inState = self.getStatePtr(self.cmodel.s,stateIds[j])
                for k in range(inState.in_states):
                    inId = ghmmwrapper.get_arrayint(inState.in_id,k)

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
        ghmmwrapper.ghmm_c_baum_welch(self.BWcontext)        
        likelihood = ghmmwrapper.get_arrayd(self.BWcontext.logp, 0)
        #(steps_made, loglikelihood_array, scale_array) = self.baumWelchStep(nrSteps,
        #                                                                    loglikelihoodCutoff)
        self.baumWelchDelete()

        return likelihood

    def baumWelchSetup(self, trainingSequences, nrSteps):
        """ Setup necessary temporary variables for Baum-Welch-reestimation.
            Use baum_welch_setup and baum_welch_step if you want more control
            over the training, compute diagnostics or do noise-insertion

            training_sequences can either be a SequenceSet or a Sequence

            Return: a C-array of type ghmm_c_baum_welch_context
        """
        self.BWcontext = ghmmwrapper.smosqd_t_array(1)
        self.BWcontext.smo  = self.cmodel
        self.BWcontext.sqd  = trainingSequences.cseq    # copy reference to ghmm_cseq
        self.BWcontext.logp = ghmmwrapper.double_array(1) # place holder for sum of log-likelihood
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

        ghmmwrapper.freearray(self.BWcontext.logp)
        self.BWcontext.logp = None
        ghmmwrapper.free_smosqd_t(self.BWcontext)
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
            state = self.getStatePtr(self.cmodel.s,i)
            pi.append(state.pi)

            B[i][0] = ghmmwrapper.get_arrayd(state.mue,0)
            B[i][1] =  ghmmwrapper.get_arrayd(state.u,0)

            for j in range(state.out_states):
                state_index = ghmmwrapper.get_arrayint(state.out_id,j)
                A[i][state_index] = ghmmwrapper.get_2d_arrayd(state.out_a,0,j)

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
        
        if not (self.cmodel.model_type & 4):   # check for kSilentStates
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
        state = ghmmwrapper.get_sstate(self.cmodel, i)
        mu = ghmmwrapper.get_arrayd(state.mue, comp)
        sigma = ghmmwrapper.get_arrayd(state.u,comp)
        weigth = ghmmwrapper.get_arrayd(state.c,comp)
        return (mu, sigma, weigth)

    def setEmission(self, i, comp,(mu, sigma, weight)):
        """ Set the emission distributionParameters for component 'comp' in state 'i'. """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = ghmmwrapper.get_sstate(self.cmodel, i)
        ghmmwrapper.set_arrayd(state.mue, comp, float(mu))  # GHMM C is german: mue instead of mu
        ghmmwrapper.set_arrayd(state.u, comp, float(sigma))
        ghmmwrapper.set_arrayd(state.c, comp, float(weight))
	
    def getEmissionProbability(self, value, state):
        return ghmmwrapper.ghmm_c_calc_b (self.cmodel, state, value)
      

    def logprob(self, emissionSequence, stateSequence):
        """ log P[ emissionSequence, stateSequence| m] """
        
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)
                        		
        state = self.getStatePtr( self.cmodel.s, stateSequence[0] )
        emissionProb = self.getEmissionProbability(emissionSequence[0],stateSequence[0])
        
    	if (emissionProb == 0): # zero ??? or some small constant?
            raise SequenceCannotBeBuild, "first symbol " + str(emissionSequence[0]) + " not emitted by state " + str(stateSequence[0])
                        
        logP = math.log(state.pi * emissionProb )
        
        #symbolIndex = 1

        try:        
            for i in range(len(emissionSequence)-1):
                cur_state = self.getStatePtr( self.cmodel.s, stateSequence[i] )
                next_state = self.getStatePtr( self.cmodel.s, stateSequence[i+1] )
		
                for j in range(cur_state.out_states):
                    out_id = ghmmwrapper.get_arrayint(cur_state.out_id,j)
                    if out_id == stateSequence[i+1]:
                        emissionProb = self.getEmissionProbability(emissionSequence[i+1],out_id)
                        #symbolIndex += 1
                        if emissionProb == 0:
                            raise SequenceCannotBeBuild, "symbol " + str(emissionSequence[i+1]) + " not emitted by state "+ str(stateSequence[i+1])
                        logP += math.log( ghmmwrapper.get_2d_arrayd(cur_state.out_a,0,j) * emissionProb)
                        break
        except IndexError:
            pass
        return logP

    def getPrior(self):
         return self.cmodel.prior

    def setPrior(self, prior):
         self.cmodel.prior = prior  
    
    def __str__(self):
        "defines string representation"
        hmm = self.cmodel

        strout = ["\nOverview of HMM:"]
        strout.append("\nNumber of states: "+ str(hmm.N))
        strout.append("\nNumber of output distributions per state: "+ str(hmm.M))

        for k in range(hmm.N):
            state = ghmmwrapper.get_sstate(hmm, k)
            strout.append("\n\nState number "+ str(k) +":")
            strout.append("\nInitial probability: " + str(state.pi))
            strout.append("\n"+ str(hmm.M) + " density function(s):\n")

            weight = ""
            mue = ""
            u =  ""

            for outp in range(hmm.M):
                weight += str(ghmmwrapper.get_arrayd(state.c,outp))+", "
                mue += str(ghmmwrapper.get_arrayd(state.mue,outp))+", "
                u += str(ghmmwrapper.get_arrayd(state.u,outp))+", "

            strout.append("  pdf component weights : " + str(weight) + "\n")
            strout.append("  mean vector: " + str(mue) + "\n")
            strout.append("  variance vector: " + str(u) + "\n")

            for c in range(self.cmodel.cos):
                strout.append("\n  Class : " + str(c)                )
                strout.append("\n    Outgoing transitions:")
                for i in range( state.out_states):
                    strout.append("\n      transition to state " + str(ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability = " + str(ghmmwrapper.get_2d_arrayd(state.out_a,c,i)))
                strout.append("\n    Ingoing transitions:")
                for i in range(state.in_states):
                    strout.append("\n      transition from state " + str(ghmmwrapper.get_arrayint(state.in_id,i) ) +" with probability = "+ str(ghmmwrapper.get_2d_arrayd(state.in_a,c,i)))


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
            state = self.getStatePtr(self.cmodel.s,i)
            pi.append(state.pi)

            B[i].append( ghmmhelper.arrayd2list(state.mue,self.cmodel.M) )
            B[i].append( ghmmhelper.arrayd2list(state.u,self.cmodel.M) )
            B[i].append( ghmmhelper.arrayd2list(state.c,self.cmodel.M) )

            for j in range(state.out_states):
                state_index = ghmmwrapper.get_arrayint(state.out_id,j)
                A[i][state_index] = ghmmwrapper.get_2d_arrayd(state.out_a,0,j)

        return [A,B,pi]


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
         return ghmmwrapper.ghmm_c_calc_b(self.cmodel, state, value)       
 
     def getEmission(self, i, comp):
        """ Return the paramenters of component 'comp' in state 'i'
        (type, mu, sigma^2, weight) - for a gaussian component
        (type, mu, sigma^2,  min, weight) - for a left trunc gaussian
        (type, mu, sigma^2, max, weight) - for a right trunc gaussian
        (type, max, mix, weight) - for a uniform
        """
        state = ghmmwrapper.get_sstate(self.cmodel, i)
        mu = ghmmwrapper.get_arrayd(state.mue, comp) 
        sigma = ghmmwrapper.get_arrayd(state.u,comp)
        weigth = ghmmwrapper.get_arrayd(state.c,comp)
        type = ghmmwrapper.get_density(state,comp)
        if ((type == ghmmwrapper.uniform) or (type == ghmmwrapper.normal)):
          return (type, mu, sigma, weigth)
        else:
          a = ghmmwrapper.get_arrayd(state.a,comp)
          return (type, mu, sigma, a, weigth)

     def setEmission(self, i, comp, type, (mu, sigma, weight, a)):
        """ Set the emission distributionParameters for component 'comp' in state 'i'. """

        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = ghmmwrapper.get_sstate(self.cmodel, i)
        ghmmwrapper.set_arrayd(state.mue, comp, float(mu))  # GHMM C is german: mue instead of mu
        ghmmwrapper.set_arrayd(state.u, comp, float(sigma))
        ghmmwrapper.set_arrayd(state.c, comp, float(weight))
        ghmmwrapper.set_arrayd(state.a, comp, float(a))        
        ghmmwrapper.set_density(state, comp, int(type))
	
     def getEmissionProbability(self, value, state):
        return ghmmwrapper.ghmm_c_calc_b (self.cmodel, state, value)
         
     def __str__(self):
         "defines string representation"
         hmm = self.cmodel
 
         strout = "\nOverview of HMM:"
         strout.append("\nNumber of states: "+ str(hmm.N))
         strout.append("\nNumber of output distributions per state: "+ str(hmm.M))
 
         for k in range(hmm.N):
             state = ghmmwrapper.get_sstate(hmm, k)
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
                 weight += str(ghmmwrapper.get_arrayd(state.c,outp))+", "
                 mue += str(ghmmwrapper.get_arrayd(state.mue,outp))+", "
                 u += str(ghmmwrapper.get_arrayd(state.u,outp))+", "
                 a += str(ghmmwrapper.get_arrayd(state.a,outp))+", "
                 density += str(ghmmwrapper.get_arrayd(state.a,outp))+","
                 c += str(ghmmwrapper.get_arrayd(state.c,outp))+","
                 
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
                     strout.append("\n      transition to state " + str(ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability = " + str(ghmmwrapper.get_2d_arrayd(state.out_a,c,i)))
                 strout.append("\n    Ingoing transitions:")
                 for i in range(state.in_states):
                     strout.append("\n      transition from state " + str(ghmmwrapper.get_arrayint(state.in_id,i) ) +" with probability = "+ str(ghmmwrapper.get_2d_arrayd(state.in_a,c,i)))
 
 
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
             state = self.getStatePtr(self.cmodel.s,i)
             pi.append(state.pi)
 
             B[i].append( ghmmhelper.arrayd2list(state.mue,self.cmodel.M) )
             B[i].append( ghmmhelper.arrayd2list(state.u,self.cmodel.M) )
             B[i].append( ghmmhelper.arrayd2list(state.c,self.cmodel.M) )
             B[i].append( ghmmhelper.arrayd2list(state.a,self.cmodel.M) )
             
             d.append([]);

             d[i].append( ghmmhelper.arrayd2list(state.density,self.cmodel.M) )
                      
 
             for j in range(state.out_states):
                 state_index = ghmmwrapper.get_arrayint(state.out_id,j)
                 A[i][state_index] = ghmmwrapper.get_2d_arrayd(state.out_a,0,j)
 
         return [A,B,pi,d]
	
def HMMDiscriminativeTraining(HMMList, SeqList, nrSteps = 50, gradient = 0):
    """ """
     
    # XXX working ? 
    
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

    HMMArray = ghmmwrapper.modelarray_alloc(inplen)
    SeqArray = ghmmwrapper.seqarray_alloc(inplen)
    
    for i in range(inplen):
        ghmmwrapper.modelarray_setptr(HMMArray, HMMList[i].cmodel, i)
        ghmmwrapper.seqarray_setptr(SeqArray, SeqList[i].cseq, i)

    ghmmwrapper.ghmm_d_discriminative(HMMArray, SeqArray, inplen, nrSteps, gradient)

    for i in range(inplen):
        HMMList[i].cmodel = ghmmwrapper.modelarray_getptr(HMMArray, i)
        SeqList[i].cseq   = ghmmwrapper.seqarray_getptr(SeqArray, i)

    ghmmwrapper.freearray(HMMArray)
    ghmmwrapper.freearray(SeqArray)

    return HMMDiscriminativePerformance(HMMList, SeqList)


        
        
def HMMDiscriminativePerformance(HMMList, SeqList):

    if len(HMMList) != len(SeqList):
        raise TypeRrror, 'Inputs not equally long'

    inplen = len(HMMList)
    
    single = [0.0] * inplen

    HMMArray = ghmmwrapper.modelarray_alloc(inplen)
    SeqArray = ghmmwrapper.seqarray_alloc(inplen)

    for i in range(inplen):
        ghmmwrapper.modelarray_setptr(HMMArray, HMMList[i].cmodel, i)
        ghmmwrapper.seqarray_setptr(SeqArray, SeqList[i].cseq, i)

    retval = ghmmwrapper.ghmm_d_discrim_performance(HMMArray, SeqArray, inplen)
    
    ghmmwrapper.freearray(HMMArray)
    ghmmwrapper.freearray(SeqArray)

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
        self.pairIndexFunction = ghmmwrapper.ghmm_dp_pair

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

class ComplexEmissionSequence:
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
                
        # necessary C functions for accessing the ghmm_dpseq struct
        self.getPtr = ghmmwrapper.get_col_pointer_int # defines C function to be used to access a single sequence
        self.getSymbol = ghmmwrapper.get_2d_arrayint
        self.setSymbol = ghmmwrapper.set_2d_arrayint
        self.cleanFunction = ghmmwrapper.ghmm_dpseq_free
        self.setDiscreteSequenceFunction = ghmmwrapper.ghmm_dpseq_set_discrete
        self.setContinuousSequenceFunction = ghmmwrapper.ghmm_dpseq_set_continuous
        self.getDiscreteSequenceFunction = ghmmwrapper.ghmm_dpseq_get_discrete
        self.getContinuousSequenceFunction = ghmmwrapper.ghmm_dpseq_get_continuous
        self.cseq = ghmmwrapper.ghmm_dpseq_init(self.length,
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
            pointerDiscrete = self.getDiscreteSequenceFunction(self.cseq, i)
            for j in range(len(self)):
                ghmmwrapper.set_arrayint(pointerDiscrete, j, internalInput[j])
            # self.setDiscreteSequenceFunction(self.cseq, i, seq)

        for i in range(len(self.continuousInputs)):
            seq = [float(x) for x in self.continuousInputs[i]]
            seq = ghmmhelper.list2arrayd(seq)
            pointerContinuous = self.getContinuousSequenceFunction(self.cseq,i)
            for j in range(len(self)):
                ghmmwrapper.set_arrayd(pointerContinuous, j, self.continuousInputs[i][j])
            # self.setContinuousSequenceFunction(self.cseq, i, seq)

    def __del__(self):
        """
        Deallocation of C sequence struct.
        """
        clean = self.cleanFunction(self.cseq, len(self.discreteDomains),
                                   len(self.continuousDomains))
        if (clean == -1):
            log.error( "could not clean sequence")
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
        int_pointer = self.getDiscreteSequenceFunction(self.cseq, index)
        internal = ghmmhelper.arrayint2list(int_pointer, len(self))
        int_pointer = None
        return internal
    
    def getInternalContinuousSequence(self, index):
        """
        access the underlying C structure and return the internal
        representation of the continuous sequence number 'index'
        @param index: number of the continuous sequence
        @return: a python list of floats
        """
        d_pointer = self.getContinuousSequenceFunction(self.cseq, index)
        internal = ghmmhelper.arrayd2list(d_pointer, len(self))
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
        
        # Assignment of the C function names to be used with this model type
        self.freeFunction = ghmmwrapper.ghmm_dp_free
        # self.samplingFunction = ghmmwrapper.ghmm_d_generate_sequences
        self.viterbiFunction = ghmmwrapper.ghmm_dp_viterbi
        self.viterbiPropagateFunction = ghmmwrapper.ghmm_dp_viterbi_propagate
        self.viterbiPropagateSegmentFunction = ghmmwrapper.ghmm_dp_viterbi_propagate_segment
        # self.forwardFunction = 
        # self.forwardAlphaFunction = 
        # self.backwardBetaFunction = 
        self.getStatePtr = ghmmwrapper.get_pstateptr 
        self.fileWriteFunction = ghmmwrapper.call_model_print
        self.getModelPtr = ghmmwrapper.get_model_ptr
        self.castModelPtr = ghmmwrapper.cast_model_ptr
        # self.distanceFunction = ghmmwrapper.ghmm_d_prob_distance
        self.logPFunction = ghmmwrapper.ghmm_dp_viterbi_logp
        self.states = {}

    def __str__(self):
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
            state = ghmmwrapper.get_pstateptr(hmm.s, k)
            strout.append("\n\nState number "+ str(k) +":")
            strout.append("\nInitial probability: " + str(state.pi))
            strout.append("\nOutput probabilites: ")
            #strout.append(str(ghmmwrapper.get_arrayd(state.b,outp)))
            strout.append("\n")

            strout.append("\nOutgoing transitions:")
            for i in range( state.out_states):
                strout.append("\ntransition to state " + str( ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.out_a,i)))
            strout.append("\nIngoing transitions:")
            for i in range(state.in_states):
                strout.append("\ntransition from state " + str( ghmmwrapper.get_arrayint(state.in_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.in_a,i)))
                strout.append("\nint fix:" + str(state.fix) + "\n")

        if hmm.model_type & 4:
            strout.append("\nSilent states: \n")
            for k in range(hmm.N):
                strout.append(str(ghmmwrapper.get_arrayint(hmm.silent,k)) + ", ")
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
        log_p_ptr = ghmmwrapper.double_array(1)
        length_ptr = ghmmwrapper.int_array(1)
        # call log_p and length will be passed by reference
        cpath = self.viterbiFunction(self.cmodel,
                                    complexEmissionSequenceX.cseq,
                                    complexEmissionSequenceY.cseq,
                                    log_p_ptr, length_ptr)
        # get the values from the pointers
        log_p = ghmmwrapper.get_arrayd(log_p_ptr, 0)
        length = ghmmwrapper.get_arrayint(length_ptr, 0)
        path = [ghmmwrapper.get_arrayint(cpath, x) for x in range(length)]
        # free the memory
        ghmmwrapper.freearray(log_p_ptr)
        ghmmwrapper.freearray(length_ptr)
        ghmmwrapper.freearray(cpath)
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
        log_p_ptr = ghmmwrapper.double_array(1)
        length_ptr = ghmmwrapper.int_array(1)
        # call log_p and length will be passed by reference
        if (not (startX and startY and stopX and stopY and startState and stopState and startLogp)):
            cpath = self.viterbiPropagateFunction(
                self.cmodel,
                complexEmissionSequenceX.cseq,
                complexEmissionSequenceY.cseq,
                log_p_ptr, length_ptr,
                self.maxSize)
        else:
            if (stopLogp == None):
                stopLogp = 0
            cpath = self.viterbiPropagateSegmentFunction(
                self.cmodel,
                complexEmissionSequenceX.cseq,
                complexEmissionSequenceY.cseq,
                log_p_ptr, length_ptr, self.maxSize,
                startX, startY, stopX, stopY, startState, stopState,
                startLogp, stopLogp)
          
        # get the values from the pointers
        log_p = ghmmwrapper.get_arrayd(log_p_ptr, 0)
        length = ghmmwrapper.get_arrayint(length_ptr, 0)
        path = [ghmmwrapper.get_arrayint(cpath, x) for x in range(length)]
        # free the memory
        ghmmwrapper.freearray(log_p_ptr)
        ghmmwrapper.freearray(length_ptr)
        ghmmwrapper.freearray(cpath)
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
        cpath = ghmmhelper.list2arrayint(path)
        logP = self.logPFunction(self.cmodel, complexEmissionSequenceX.cseq,
                                 complexEmissionSequenceY.cseq,
                                 cpath, len(path))
        ghmmwrapper.freearray(cpath)
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
        self.cmodel.size_of_alphabet = ghmmhelper.list2arrayint(self.alphabetSizes)

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
                c_state = ghmmwrapper.get_pstateptr(self.cmodel.s,
                                                    orders[state.index])
                for out in range(c_state.out_states):
                    outSum += ghmmwrapper.get_2d_arrayd(c_state.out_a,
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
        sizes = [alphabets[alphabet].size() for alphabet in alphabets.keys()]
        cmodel.size_of_alphabet = ghmmhelper.list2arrayint(sizes)

        # set number of d_seqs to zero. If you want to use them you have to
        # set them manually
        cmodel.number_of_d_seqs = 0

        # c array of states allocated
        cstates = ghmmwrapper.arraypstate(cmodel.N)
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
            cstate = ghmmwrapper.get_pstateptr(cstates,i)
            pystate = pystates[i]
            size = pystate.itsHMM.hmmAlphabets[pystate.alphabet_id].size()
            if (pystate.offsetX != 0 and pystate.offsetY != 0):
                size = size**2
            if (len(B[i]) != size):
                raise InvalidModelParameters("in state %s len(emissions) = %i size should be %i" % (pystate.id, len(B[i]), size))
            cstate.b = ghmmhelper.list2arrayd(B[i])
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
                cstate.out_id = ghmmhelper.list2arrayint(myoutid)
                (cstate.out_a, col_len) = ghmmhelper.list2matrixd(outprobs)
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
                cstate.in_id = ghmmhelper.list2arrayint(myinid)
                (cstate.in_a, col_len) = ghmmhelper.list2matrixd(inprobs)
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
        cmodel.silent = ghmmhelper.list2arrayint(silent_states)        
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
