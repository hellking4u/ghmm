#!/usr/bin/env python2.3F
################################################################################
#
#       This file is part of the General Hidden Markov Model Library,
#       GHMM version __VERSION__, see http://ghmm.org
#
#       file:    ghmm.py
#       authors: Benjamin Georgi, Wasinee Rungsarityotin, Alexander Schliep
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

from os import path
from math import log,ceil


print "*** I'm the ghmm in /amd/rindt/1/home/abt_vin/georgi/hmm/ghmm/ghmmwrapper ***"

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

    def __str__(self):
        strout = "GHMM Alphabet:\n"
        strout += "Number of symbols: " + str(len(self)) + "\n"
        strout += "External: " + str(self.listOfCharacters) + "\n"
        strout += "Internal: " + str(range(len(self))) + "\n"
        return strout
    

    def __len__(self):
        return len(self.listOfCharacters)
        

    def size(self):
        """ Deprecated """
        print "Warning: The use of .size() is deprecated. Use len() instead."
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

            Raises KeyError
        """
        if internal < 0 or len(self.listOfCharacters) <= internal:
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


class ContinousDistribution(Distribution):
    pass

class GaussianDistribution(ContinousDistribution):
    # XXX attributes unused at this point
    def __init__(self, domain):
        self.emissionDomain = domain
        self.mu = None
        self.sigma = None

    def set(self, (mu, sigma)):
        self.mu = mu
        self.sigma = sigma

    def get(self):
        return (self.mu, self.sigma2)


class MixtureContinousDistribution(ContinousDistribution):
    pass

class GaussianMixtureDistribution(MixtureContinousDistribution):
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


#-------------------------------------------------------------------------------
#Sequence, SequenceSet and derived  ------------------------------------------

class EmissionSequence:
    """ An EmissionSequence contains the *internal* representation of
        a sequence of emissions. It also contains a reference to the
        domain where the emission orginated from.
    """

    def __init__(self, emissionDomain, sequenceInput, labelDomain = None, labelInput = None):

        self.emissionDomain = emissionDomain

        if self.emissionDomain.CDataType == "int": # underlying C data type is integer

            # necessary C functions for accessing the sequence_t struct
            self.getPtr = ghmmwrapper.get_col_pointer_int # defines C function to be used to access a single sequence
            self.getSymbol = ghmmwrapper.get_2d_arrayint
            self.setSymbol = ghmmwrapper.set_2d_arrayint
            self.freeFunction = ghmmwrapper.call_sequence_free
            self.addSeqFunction = ghmmwrapper.sequence_add

            #create a sequence_t with state_labels, if the appropiate parameters are set
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

                self.cseq = ghmmwrapper.sequence_calloc(1)
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
                self.cseq = ghmmwrapper.sequence_calloc(1)
                self.cseq.seq = seq
                self.cseq.seq_number = 1

                # deactivating labels
                self.cseq.state_labels = None
                self.cseq.state_labels_len = None                
                
                ghmmwrapper.set_arrayint(self.cseq.seq_len,0,l[0]) 

        
            elif (isinstance(sequenceInput, str) or isinstance(sequenceInput, unicode)): # from file

                # reads in the first sequence struct in the input file
                if  not path.exists(sequenceInput):
                     raise IOError, 'File ' + str(sequenceInput) + ' not found.'
                else:
                    self.cseq  = ghmmwrapper.seq_read(sequenceInput)

            elif isinstance(sequenceInput, ghmmwrapper.sequence_t):# internal use
                if sequenceInput.seq_number > 1:
                    raise badCPointer, "Use SequenceSet for multiple sequences."
                self.cseq = sequenceInput
                if labelDomain != None:
                    self.labelDomain = labelDomain
                
                
            else:
                raise UnknownInputType, "inputType " + str(type(sequenceInput)) + " not recognized."

        elif self.emissionDomain.CDataType == "double": # underlying C data type is double

            # necessary C functions for accessing the sequence_d_t struct
            self.getPtr = ghmmwrapper.get_col_pointer_d # defines C function to be used to access a single sequence
            self.getSymbol = ghmmwrapper.get_2d_arrayd
            self.setSymbol = ghmmwrapper.set_2d_arrayd
            self.freeFunction = ghmmwrapper.call_sequence_d_free
            self.addSeqFunction = ghmmwrapper.sequence_d_add

            if isinstance(sequenceInput, list):
                (seq,l) = ghmmhelper. list2matrixd([sequenceInput])
                self.cseq = ghmmwrapper.sequence_d_calloc(1)
                self.cseq.seq = seq
                self.cseq.seq_number = 1
                ghmmwrapper.set_arrayint(self.cseq.seq_len,0,l[0])

            elif isinstance(sequenceInput, str) or isinstance(sequenceInput, unicode): # from file
                # reads in the first sequence struct in the input file
                if  not path.exists(sequenceInput):
                     raise IOError, 'File ' + str(sequenceInput) + ' not found.'
                else:
                    self.cseq  = ghmmwrapper.seq_d_read(sequenceSetInput)


            elif isinstance(sequenceInput, ghmmwrapper.sequence_d_t): # internal use
                if sequenceInput.seq_number > 1:
                    raise badCPointer, "Use SequenceSet for multiple sequences."
                self.cseq = sequenceInput
            else:
                raise UnknownInputType, "inputType " + str(type(sequenceInput)) + " not recognized."


        else:
            raise NoValidCDataType, "C data type " + str(self.emissionDomain.CDataType) + " invalid."

    def __del__(self):
        "Deallocation of C sequence struct."
        print "__del__ EmissionSequence " + str(self.cseq)
        self.freeFunction( self.cseq)
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
        if (self.emissionDomain.CDataType != "double"):
            print "WARNING: Discrete sequences do not support sequence labels."
        else:
            return ghmmwrapper.get_arrayl(self.cseq.seq_label,0)

    def setSeqLabel(self,value):
        if (self.emissionDomain.CDataType != "double"):
            print "WARNING: Discrete sequences do not support sequence labels."
        else:     
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
        strout = ""
        strout += "\nEmissionSequence Instance:\nlength " + str(ghmmwrapper.get_arrayint(seq.seq_len,0))+ ", weight " + str(ghmmwrapper.get_arrayd(seq.seq_w,0))  + ":\n"
        for j in range(ghmmwrapper.get_arrayint(seq.seq_len,0) ):
            strout += str( self.emissionDomain.external(self[j]) )   
            if self.emissionDomain.CDataType == "double":
                strout += " "

        # checking for labels
        if self.emissionDomain.CDataType == "int" and self.cseq.state_labels != None:            
            strout += "\nState labels:\n"
            for j in range(ghmmwrapper.get_arrayint(seq.state_labels_len,0) ):
                strout += str( self.labelDomain.external(ghmmwrapper.get_2d_arrayint(seq.state_labels,0,j)))+ ", "

    	return strout

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
            ghmmwrapper.call_sequence_print(fileName, self.cseq)
        if self.emissionDomain.CDataType == "double":
            ghmmwrapper.call_sequence_d_print(fileName, self.cseq,0)

    
    def setWeight(self, value):
        ghmmwrapper.set_arrayd(self.cseq.seq_w,0,value)
        self.cseq.total_w = value
        
    def getWeight(self):
        return ghmmwrapper.get_arrayd(self.cseq.seq_w,0)


class SequenceSet:
    def __init__(self, emissionDomain, sequenceSetInput, labelDomain = None, labelInput = None):
        self.emissionDomain = emissionDomain
        self.cseq = None
        
        if self.emissionDomain.CDataType == "int": # underlying C data type is integer
            # necessary C functions for accessing the sequence_t struct
            self.sequenceAllocationFunction = ghmmwrapper.sequence_calloc
            self.getPtr = ghmmwrapper.get_col_pointer_int # defines C function to be used to access a single sequence
            self.castPtr = ghmmwrapper.cast_ptr_int # cast int* to int** pointer
            self.emptySeq = ghmmwrapper.int_2d_array_nocols # allocate an empty int**
            self.setSeq = ghmmwrapper.set_2d_arrayint_col # assign an int* to a position within an int**
            self.freeFunction = ghmmwrapper.call_sequence_free
            self.addSeqFunction = ghmmwrapper.sequence_add # add sequences to the underlying C struct
            self.getSymbol = ghmmwrapper.get_2d_arrayint
            self.setSymbolSingle = ghmmwrapper.set_arrayint
            
            if (isinstance(sequenceSetInput, str)  and labelInput == None): # from file
                # reads in the first sequence struct in the input file
                if  not path.exists(sequenceSetInput):
                     raise IOError, 'File ' + str(sequenceSetInput) + ' not found.'
                else:
                     self.cseq  = ghmmwrapper.seq_read(sequenceSetInput)

            #generate a a labeled sequenceSet from a list of lists (sequences), emissionDomain, a list of list (labels) and a model
            elif isinstance(sequenceSetInput, list) and isinstance(labelInput, list) and isinstance(labelDomain, LabelDomain): 
                assert len(sequenceSetInput)==len(labelInput)

                self.labelDomain = labelDomain                
                seq_nr = len(sequenceSetInput)
                self.cseq = ghmmwrapper.sequence_calloc(seq_nr)
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
                self.cseq = ghmmwrapper.sequence_calloc(seq_nr)
                self.cseq.seq_number = seq_nr
                
                internalInput = []
                for i in range(seq_nr):
                    internalInput.append( map( self.emissionDomain.internal, sequenceSetInput[i]))

                (seq,lenghts) = ghmmhelper.list2matrixint(internalInput)
                
                self.cseq.seq = seq
                for i in range(seq_nr):
                    ghmmwrapper.set_arrayint(self.cseq.seq_len,i ,lenghts[i])


            elif isinstance(sequenceSetInput, ghmmwrapper.sequence_t): # inputType == sequence_t*
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
            # necessary C functions for accessing the sequence_d_t struct
            self.sequenceAllocationFunction =  ghmmwrapper.sequence_d_calloc
            self.getPtr = ghmmwrapper.get_col_pointer_d # defines C function to be used to access a single sequence
            self.castPtr = ghmmwrapper.cast_ptr_d # cast double* to double** pointer
            self.emptySeq = ghmmwrapper.double_2d_array_nocols  # cast double* to int** pointer
            self.setSeq = ghmmwrapper.set_2d_arrayd_col # assign a double* to a position within a double**
            self.freeFunction = ghmmwrapper.call_sequence_d_free
            self.addSeqFunction = ghmmwrapper.sequence_d_add # add sequences to the underlying C struct
            self.getSymbol = ghmmwrapper.get_2d_arrayd
            self.setSymbolSingle = ghmmwrapper.set_arrayd

                        
            if isinstance(sequenceSetInput, list): 
                seq_nr = len(sequenceSetInput)
                self.cseq = ghmmwrapper.sequence_d_calloc(seq_nr)
                self.cseq.seq_number = seq_nr

                (seq,lenghts) = ghmmhelper.list2matrixd(sequenceSetInput)
                self.cseq.seq = seq
                for i in range(seq_nr):
                    ghmmwrapper.set_arrayint(self.cseq.seq_len, i, lenghts[i])

            elif isinstance(sequenceSetInput, str) or isinstance(sequenceSetInput, unicode): # from file
                print "fromFile", sequenceSetInput
                if  not path.exists(sequenceSetInput):
                     raise IOError, 'File ' + str(sequenceSetInput) + ' not found.'
                else:
                    #self.cseq  = ghmmwrapper.seq_d_read(sequenceSetInput)
                    i = ghmmwrapper.int_array(1)
                    self.__s = ghmmwrapper.sequence_d_read(sequenceSetInput,i)
                    self.cseq = ghmmwrapper.get_seq_d_ptr(self.__s,0)

                                                     
            elif isinstance(sequenceSetInput, ghmmwrapper.sequence_d_t): # i# inputType == sequence_d_t**, internal use

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
        
        print "** SequenceSet.__del__ " + str(self.cseq)
        self.freeFunction(self.cseq)
        self.cseq = None
    
    def __str__(self):
        "Defines string representation."
        seq = self.cseq
        strout =  "\nNumber of sequences: " + str(seq.seq_number)
        
        if self.emissionDomain.CDataType == "int":
           for i in range(seq.seq_number):
                strout += "\nSeq " + str(i)+ ", length " + str(ghmmwrapper.get_arrayint(seq.seq_len,i))+ ", weight " + str(ghmmwrapper.get_arrayd(seq.seq_w,i))  + ":\n"
                for j in range(ghmmwrapper.get_arrayint(seq.seq_len,i) ):
                    strout += str( self.emissionDomain.external(( ghmmwrapper.get_2d_arrayint(self.cseq.seq, i, j) )) )
                    if self.emissionDomain.CDataType == "double":
                       strout += " "

                # checking for labels 
                if self.emissionDomain.CDataType == "int" and self.cseq.state_labels != None:            
                    strout += "\nState labels:\n"
                    for j in range(ghmmwrapper.get_arrayint(seq.state_labels_len,i) ):
                        strout += str( self.labelDomain.external(ghmmwrapper.get_2d_arrayint(seq.state_labels,i,j))) +", "

        if self.emissionDomain.CDataType == "double":
            for i in range(seq.seq_number):
                strout += "\nSeq " + str(i)+ ", length " + str(ghmmwrapper.get_arrayint(seq.seq_len,i))+ ", weight " + str(ghmmwrapper.get_arrayd(seq.seq_w,i))  + ":\n"
                for j in range(ghmmwrapper.get_arrayint(seq.seq_len,i) ):
                    strout += str( self.emissionDomain.external(( ghmmwrapper.get_2d_arrayd(self.cseq.seq, i, j) )) ) + " "

        return strout
    


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
        """ Return the weight of sequence i. Weights are used in Baum-Welch"""
        return ghmmwrapper.set_arrayd(self.cseq.seq_w,i,w)
        
    def __getitem__(self, index):
        """ Return an EmissionSequence object initialized with a reference to 
        sequence 'index'.
        
        """
        seq = self.sequenceAllocationFunction(1)
        seq.seq = self.castPtr(self.__array[index]) # int* -> int** reference
        ghmmwrapper.set_arrayint(seq.seq_len,0,ghmmwrapper.get_arrayint(self.cseq.seq_len,index))
        seq.seq_number = 1

        return EmissionSequence(self.emissionDomain, seq)        

    def getSeqLabel(self,index):
        if (self.emissionDomain.CDataType != "double"):
            print "WARNING: Discrete sequences do not support sequence labels."
        else:
            return ghmmwrapper.get_arrayl(self.cseq.seq_label,index)

    def setSeqLabel(self,index,value):
        if (self.emissionDomain.CDataType != "double"):
            print "WARNING: Discrete sequences do not support sequence labels."
        else:     
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
        seq.seq = self.emptySeq(seqNumber)
        
        # checking for state labels in the source C sequence struct
        if self.emissionDomain.CDataType == "int" and self.cseq.state_labels is not None:
            
            print "found labels !"
            seq.state_labels = ghmmwrapper.int_2d_array_nocols(seqNumber)
            seq.state_labels_len = ghmmwrapper.int_array(seqNumber)

        for i in range(seqNumber):
            
            self.setSeq(seq.seq,i,self.__array[ seqIndixes[i] ])
            ghmmwrapper.set_arrayint(seq.seq_len,i,ghmmwrapper.get_arrayint(self.cseq.seq_len,seqIndixes[i]))
            
            # Above doesnt copy seq_id or seq_label or seq_w
            seq_id = int(ghmmwrapper.get_arrayd(self.cseq.seq_id, seqIndixes[i]))
            ghmmwrapper.set_arrayd(seq.seq_id, i, seq_id)
            seq_label = ghmmwrapper.get_arrayl(self.cseq.seq_label, i)
            ghmmwrapper.set_arrayl(seq.seq_label, i, int(seq_label))

            seq_w = ghmmwrapper.get_arrayd(self.cseq.seq_w, i)
            ghmmwrapper.set_arrayd(seq.seq_w, i, seq_w)
             
            # setting labels if appropriate
            if self.emissionDomain.CDataType == "int" and self.cseq.state_labels is not None:
                self.setSeq(seq.state_labels,i, ghmmwrapper.get_col_pointer_int( self.cseq.state_labels,seqIndixes[i] ) )
                ghmmwrapper.set_arrayint(seq.state_labels_len,i,ghmmwrapper.get_arrayint(self.cseq.state_labels_len, seqIndixes[i]) )            

        seq.seq_number = seqNumber
        
        return SequenceSet(self.emissionDomain, seq)
        
    def write(self,fileName):
        "Writes (appends) the SequenceSet into file 'fileName'."
        
        # different function signatures require explicit check for C data type
        if self.emissionDomain.CDataType == "int":
            ghmmwrapper.call_sequence_print(fileName, self.cseq)
        if self.emissionDomain.CDataType == "double":    
            ghmmwrapper.call_sequence_d_print(fileName, self.cseq,0)


def SequenceSetOpen(emissionDomain, fileName):
    """ Reads a sequence file with multiple sequence sets. 

    Returns a list of SequenceSet objects.
    
    """

    if not path.exists(fileName):
        raise IOError, 'File ' + str(fileName) + ' not found.'

    
    if emissionDomain.CDataType == "int":
        readFile = ghmmwrapper.sequence_read
        seqPtr = ghmmwrapper.get_seq_ptr
    elif emissionDomain.CDataType == "double":
        readFile = ghmmwrapper.sequence_d_read
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
            if not path.exists(fileName):
                raise IOError, 'File ' + str(fileName) + ' not found.'

    	if self.defaultFileType == GHMM_FILETYPE_XML: # XML file
    	    hmm_dom = xmlutil.HMM(fileName)
    	    emission_domain = hmm_dom.AlphabetType()
    	    if emission_domain == int:
                [alphabets, A, B, pi, state_orders] = hmm_dom.buildMatrices()

    		emission_domain = Alphabet(alphabets)
    		distribution = DiscreteDistribution(emission_domain)
    		# build adjacency list

                # check for background distributions
                (background_dist, orders, code2name) = hmm_dom.getBackgroundDist()
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
                     print "DEBUG model_type %x" % m.cmodel.model_type
                     print "DEBUG background_id", ghmmhelper.arrayint2list(m.cmodel.background_id, m.N)
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
                    print "DEBUG durations: ", durations
                    m.extendDurations(durations)
                        
	    	return m
	    
        elif self.defaultFileType == GHMM_FILETYPE_SMO:
	        # MO & SMO Files, format is deprecated
    	    (hmmClass, emission_domain, distribution) = self.determineHMMClass(fileName)

            #print "determineHMMClass = ",  (hmmClass, emission_domain, distribution)

            nrModelPtr = ghmmwrapper.int_array(1)
    	    
            # XXX broken since silent states are not supported by .smo file format 
            if hmmClass == DiscreteEmissionHMM:
                print nrModelPtr
                models = ghmmwrapper.model_read(fileName, nrModelPtr)
                getPtr = ghmmwrapper.get_model_ptr
            else:
                models = ghmmwrapper.smodel_read(fileName, nrModelPtr)
                getPtr = ghmmwrapper.get_smodel_ptr

            nrModels = ghmmwrapper.get_arrayint(nrModelPtr, 0)
            print nrModels
            if modelIndex == None:
                cmodel = getPtr(models, 0) 
                print cmodel
            else:
                if modelIndex < nrModels:
                    cmodel = getPtr(models, modelIndex) 
                else:
		            return None

            m = hmmClass(emission_domain, distribution(emission_domain), cmodel)
            return m
        
        else:   
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

    
    def all(self, fileName):
        
        if not path.exists(fileName):
          raise IOError, 'File ' + str(fileName) + ' not found.'

        # MO & SMO Files
        (hmmClass, emission_domain, distribution) = self.determineHMMClass(fileName)
        nrModelPtr = ghmmwrapper.int_array(1)
        if hmmClass == DiscreteEmissionHMM:
            models = hmmwrapper.model_read(fileName, nrModelPtr)
            getPtr = ghmmwrapper.get_model_ptr
        else:
            models = ghmmwrapper.smodel_read(fileName, nrModelPtr)
            getPtr = ghmmwrapper.get_smodel_ptr

        nrModels = ghmmwrapper.get_arrayint(nrModelPtr, 0)
        result = []
        for i in range(nrModels):
            cmodel = getPtr(models, i)
            m = hmmClass(emission_domain, distribution(emission_domain), cmodel)
            result.append(m)
        return result
        

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
                        print "HMMOpenFactory:determineHMMClass: both HMM and SHMM?", emission_domain
                    else:
                        emission_domain = 'int'
                    
                if shmm != None:
                    if emission_domain != None and emission_domain != 'double':
                        print "HMMOpenFactory:determineHMMClass: both HMM and SHMM?", emission_domain
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
            if cos == 1 and M == 1 and density == 0:
                emission_domain = Float()
                distribution = GaussianDistribution
                hmm_class = GaussianEmissionHMM            
                return (hmm_class, emission_domain, distribution)

            elif cos == 1 and M > 1 and density == 0:
                emission_domain = Float()
                distribution = GaussianMixtureDistribution
                hmm_class = GaussianMixtureHMM
                return (hmm_class, emission_domain, distribution)

        return (None, None, None)
            
HMMOpenHMMER = HMMOpenFactory(GHMM_FILETYPE_HMMER) # read single HMMER model from file
HMMOpen = HMMOpenFactory(GHMM_FILETYPE_SMO)
HMMOpenXML = HMMOpenFactory(GHMM_FILETYPE_XML)


def readMultipleHMMERModels(fileName):
    """
        Reads a file containing multiple HMMs in HMMER format, returns list of
        HMM objects.

    """
    
    if not path.exists(fileName):
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
            print "reading model " + name + " "
            
        match = res.match(line)
        if match:
            fileLike = StringIO.StringIO(string)
            modelList.append(HMMOpenHMMER(fileLike))
            string = ""
            match = None

    return modelList
    

class HMMFromMatricesFactory(HMMFactory):
    def __call__(self, emissionDomain, distribution, A, B, pi, hmmName = None, labelDomain= None, labelList = None):
        if isinstance(emissionDomain,Alphabet):
            
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

                maxorder = 0

                #initialize states
                for i in range(cmodel.N):
                    state = ghmmwrapper.get_stateptr(states,i)
                    # compute state order
                    if cmodel.M > 1:
                        order = ( log(len(B[i]),cmodel.M) -1)
                    else:
                        order = len(B[i]) - 1
                        
                    #print "order in state ",i," = ", order
                    # check or valid number of emission parameters
                    order = int(order)
                    if  cmodel.M**(order+1) == len(B[i]):
                        state.order = order
                    else:
                        raise InvalidModelParameters, "The number of "+str(len(B[i]))+ " emission parameters for state "+str(i)+" is invalid. State order can not be determined."
                    
                    
                    state.b = ghmmhelper.list2arrayd(B[i])
                    
                    
                    state.pi = pi[i]
                    if state.order > maxorder:
                        print "state ",i,", order ",state.order
                        
                        maxorder = state.order
                    
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
                    #fix probabilities by reestimation, else 0
                    state.fix = 0

                cmodel.s = states
                cmodel.model_type += silent_flag
                cmodel.silent = ghmmhelper.list2arrayint(silent_states)

                cmodel.maxorder = maxorder
                if cmodel.maxorder > 0:
                    #print "Set kHigherOrderEmissions"
                    cmodel.model_type += 16     #kHigherOrderEmissions

                # initialize lookup table for powers of the alphabet size,
                # speeds up models with higher order states
                powLookUp = [1] * (maxorder+2)
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
                    ghmmwrapper.smodel_class_change_alloc(cmodel)
                else: 
                    cos = 1
                    cmodel.N = len(A)
                    A = [A]
                
                cmodel.cos = ghmmhelper.classNumber(A)  # number of transition classes in GHMM

                print "cmodel.cos = ",cmodel.cos

                states = ghmmwrapper.arraysstate(cmodel.N)

                # Switching functions and transition classes are handled
                # elswhere

                #initialize states
                for i in range(cmodel.N):

                    #print "  state " + str(i) + ":"
                    
                    state = ghmmwrapper.get_sstate_ptr(states,i)
                    state.pi = pi[i]

                    # allocate arrays of emmission parameters
                    state.c = ghmmhelper.list2arrayd([1.0]) # Mixture weights. Unused
                    (mu, sigma) = B[i]
                    state.mue = ghmmhelper.list2arrayd([mu]) #mu = mue in GHMM C-lib.
                    state.u = ghmmhelper.list2arrayd([sigma])
                    
                    # mixture fixing deactivated by default
                    state.mixture_fix = ghmmhelper.list2arrayint([0])
                    
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
                print " ** mixture model"
                
                # Interpretation of B matrix for the mixture case (Example with three states and two components each):
                #  B = [ 
                #      [ ["mu11","mu12"],["sig11","sig12"],["w11","w12"]   ],
                #      [  ["mu21","mu22"],["sig21","sig22"],["w21","w22"]  ],
                #      [  ["mu31","mu32"],["sig31","sig32"],["w31","w32"]  ],
                #      ]
                
                cmodel = ghmmwrapper.new_smodel()
                cmodel.N = len(A)
                cmodel.M = len(B[0][0]) # Number of mixture componenents for emission distribution
                cmodel.prior = -1 # Unused
                
                # determining number of transition classes
                cos = 0
                if type(A[0][0]) == list:
                    cos = len(A)
                else: 
                    cos = 1
                    A = [A]
                cmodel.cos = ghmmhelper.classNumber(A)  # number of transition classes in GHMM
                
                states = ghmmwrapper.arraysstate(cmodel.N)

                # Switching functions and transition classes are handled
                # elswhere

                #initialize states
                for i in range(cmodel.N):

                    print "  state " + str(i) + ":"
                    
                    state = ghmmwrapper.get_sstate_ptr(states,i)
                    state.pi = pi[i]

                    # allocate arrays of emmission parameters
                    state.c = ghmmhelper.list2arrayd([1.0]) # Mixture weights. Unused
                    mu_list = B[i][0]
                    sigma_list = B[i][1]
                    weight_list = B[i][2]
                    
                    state.mue = ghmmhelper.list2arrayd(mu_list) #mu = mue in GHMM C-lib.
                    state.u = ghmmhelper.list2arrayd(sigma_list)
                    state.c = ghmmhelper.list2arrayd(weight_list)

                    # mixture fixing deactivated by default
                    state.mixture_fix = ghmmhelper.list2arrayint([0] * cmodel.M)
                    
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
                

            else:
                raise GHMMError(type(distribution),
                                "Cannot construct model for this domain/distribution combination") 


HMMFromMatrices = HMMFromMatricesFactory()

#-------------------------------------------------------------------------------
#- Background distribution

class BackgroundDistribution:
    def __init__(self,emissionDomain, bgInput):
        
        if type(bgInput) == list:
            self.emissionDomain = emissionDomain
            distNum = len(bgInput)
        
            order = ghmmwrapper.int_array(distNum)
            b = ghmmwrapper.double_2d_array_nocols(distNum)
            for i in range(distNum):
                if len(emissionDomain) > 1:
                    o = log(len(bgInput[i]),len(emissionDomain)) -1
                else:
                    o = len(bgInput[i]) - 1
                         
                assert (o % 1) == 0, "Ivalid order of distribution " + str(i) + ": " + str(o)

                ghmmwrapper.set_arrayint(order,i, int(o))
                # dynamic allocation, rows have different lenghts
                b_i = ghmmhelper.list2arrayd(bgInput[i])
                ghmmwrapper.set_2d_arrayd_col(b,i,b_i)
    
            self.cbackground = ghmmwrapper.model_alloc_background_distributions(distNum,len(emissionDomain) ,order, b)

        elif isinstance(bgInput,ghmmwrapper.background_distributions ):
            self.cbackground = bgInput
            self.emissionDomain = emissionDomain
            
        else:
            raise TypeError, "Input type "+str(type(bgInput)) +" not recognized."    

    def __del__(self):
        print "** Freeing ", self.cbackground
        ghmmwrapper.model_free_background_distributions(self.cbackground)
        self.cbackground = None
    
    def __str__(self):
        outstr = "BackgroundDistribution instance:\n"
        outstr += "Number of distributions: " + str(self.cbackground.n)+"\n\n"
        outstr += str(self.emissionDomain) + "\n"
        d = ghmmhelper.matrixd2list(self.cbackground.b,self.cbackground.n,len(self.emissionDomain))
        outstr += "Distributions:\n"   
        for i in range(self.cbackground.n):
            outstr += "  Order: " + str(ghmmwrapper.get_arrayint(self.cbackground.order,i))+"\n"
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

        self.silent = 0   # flag for silent states
        
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
        #print "HMM.__del__", self.cmodel
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
                
                print "Warning: forward returned -1: Sequence", i,"cannot be build."
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

        ghmmwrapper.free_arrayd(likelihood)
        likelihood = None
        return likelihoodList

    ## Further Marginals ...

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
        
    # extern double smodel_prob_distance(smodel *cm0, smodel *cm, int maxT, int symmetric, int verbose);
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
        calpha = ghmmwrapper.matrix_d_alloc(t,self.N)
        cscale = ghmmwrapper.double_array(t)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)

        error = self.forwardAlphaFunction(self.cmodel, seq,t, calpha, cscale, unused)
        if error == -1:
            print "ERROR: forward finished with -1: EmissionSequence cannot be build."
            
        
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
            print "ERROR: backward finished with -1: EmissionSequence cannot be build."
            
        
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

            if seq_len != 0:
                viterbiPath = self.viterbiFunction(self.cmodel,seq,seq_len,log_p)
            else:
                viterbiPath = None

            onePath = []
            
            # for model types without possible silent states the length of the viterbi path is known
            if self.silent == 0:            
                for j in range(seq_len):                
                    onePath.append(ghmmwrapper.get_arrayint(viterbiPath,j))
            
            # in the silent case we have to reversely append as long as the path is positive because unused positions
            # are initialised with -1 on the C level.
            elif self.silent == 1:
                
                for j in range( ( seq_len * self.N )-1,-1,-1): # maximum length of a viterbi path for a silent model
                    d = ghmmwrapper.get_arrayint(viterbiPath,j)
                                   
                    if d >= 0:
                        onePath.insert(0,d)
                    else:
                        break

            allPaths.append(onePath)
            allLogs.append(ghmmwrapper.get_arrayd(log_p, 0))
        
        ghmmwrapper.free_arrayd(log_p)
        ghmmwrapper.free_arrayi(viterbiPath)  
        viterbiPath = None
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

        ghmmwrapper.model_set_transition(self.cmodel, i, j, prob)


    def getEmission(self, i):
        """ Accessor function for the emission distribution parameters of state 'i'.

            For discrete models the distribution over the symbols is returned,
            for continous models a matrix of the form
            [ [mu_1, sigma_1, weight_1] ... [mu_M, sigma_M, weight_M]  ] is returned.

        """
        if self.emissionDomain.CDataType == "int": # discrete emissions.
            state = self.getStatePtr(self.cmodel.s,i)
            emissions = ghmmhelper.arrayd2list(state.b,self.M**(state.order+1))
            return emissions

        elif self.emissionDomain.CDataType == "double": # continous emissions
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
        strout = ""
        first = 0
        if model_type &  2:         #kLeftRight
            strout += "kLeftRight "
        if model_type &  4:         #kSilentStates
            strout += "kSilentStates "
        if model_type &  8:         #kTiedEmissions
            strout += "kTiedEmissions "
        if model_type & 16:         #kHigherOrderEmissions
            strout += "kHigherOrderEmissions "
        if model_type & 32:         #kHasBackgroundDistributions
            strout += "kHasBackgroundDistributions "
        if model_type & 64:         #kClassLabels
            strout += "kClassLabels "
        if model_type == 0:         #kNotSpecified
            strout = "kNotSpecified"
        return strout


def HMMwriteList(fileName,hmmList):
    if path.exists(fileName):
        print "HMMwriteList warning: File " + str(fileName) + " already exists. " + str(len(hmmList)) + " new models will be appended."
    for model in hmmList:
        model.write(fileName)

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
        self.silent = 1  # flag indicating whether the model type does include silent states

        self.model_type = self.cmodel.model_type  # model type
        self.maxorder = self.cmodel.maxorder
        self.background = None
        
        # Assignment of the C function names to be used with this model type
        self.freeFunction = ghmmwrapper.call_model_free
        self.samplingFunction = ghmmwrapper.model_generate_sequences
        self.viterbiFunction = ghmmwrapper.viterbi
        self.forwardFunction = ghmmwrapper.foba_forward_lean
        self.forwardAlphaFunction = ghmmwrapper.foba_forward      
        self.backwardBetaFunction = ghmmwrapper.foba_backward
        self.getStatePtr = ghmmwrapper.get_stateptr 
        self.fileWriteFunction = ghmmwrapper.call_model_print
        self.getModelPtr = ghmmwrapper.get_model_ptr
        self.castModelPtr = ghmmwrapper.cast_model_ptr
        self.distanceFunction = ghmmwrapper.model_prob_distance
    
    def __del__(self):
        if self.cmodel.tied_to is not None:
            self.removeTiegroups()
        HMM.__del__(self)
        
    def __str__(self):
        hmm = self.cmodel
        strout = "\nGHMM Model\n"
        strout += "Name: " + str(self.cmodel.name)
        strout += "\nModelflags: "+ self.printtypes(self.cmodel.model_type)
        strout += "\nNumber of states: "+ str(hmm.N)
        strout += "\nSize of Alphabet: "+ str(hmm.M)
        for k in range(hmm.N):
            state = ghmmwrapper.get_stateptr(hmm.s, k)
            strout += "\n\nState number "+ str(k) +":"
            strout += "\nState order: " + str(state.order)
            strout += "\nInitial probability: " + str(state.pi)
            #strout += "\nsilent state: " + str(get_arrayint(self.model.silent,k))
            strout += "\nOutput probabilites: "
            for outp in range(hmm.M**(state.order+1)):
                strout+=str(ghmmwrapper.get_arrayd(state.b,outp))
                if outp % hmm.M == hmm.M-1:
                    strout += "\n"
                else:
                    strout += ", "

            strout += "\nOutgoing transitions:"
            for i in range( state.out_states):
                strout += "\ntransition to state " + str( ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.out_a,i))
            strout +=  "\nIngoing transitions:"
            for i in range(state.in_states):
                strout +=  "\ntransition from state " + str( ghmmwrapper.get_arrayint(state.in_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.in_a,i))
                strout += "\nint fix:" + str(state.fix) + "\n"
        strout += "\nSilent states: \n"
        for k in range(hmm.N):
            strout += str(ghmmwrapper.get_arrayint(self.cmodel.silent,k)) + ", "
        strout += "\n"
        return strout
    

    def extendDurations(self, durationlist):
        """ extend states with durations larger one
            this done by explicit state copying in C """

        for i in range(len(durationlist)):
            if durationlist[i] > 1:
                error = ghmmwrapper.model_apply_duration(self.cmodel, i, durationlist[i])
                self.N = self.cmodel.N
                if error:
                    print "ERROR: durations not applied"

    def setEmission(self, i, distributionParameters):
        """ Set the emission distribution parameters for a discrete model."""
        assert len(distributionParameters) == self.M
        # ensure proper indices
        assert 0 <= i < self.N, "Index " + str(i) + " out of bounds."

        state = self.getStatePtr(self.cmodel.s,i)

        # updating silent flag if necessary        
        if sum(distributionParameters) > 0:
            silentFlag =  0 
        else:
            silentFlag =  1 
        oldFlag = ghmmwrapper.get_arrayint(self.cmodel.silent,i)
        
        if silentFlag != oldFlag:
            ghmmwrapper.set_arrayint(self.cmodel.silent,i,silentFlag)            
            # checking if model type changes
            if silentFlag == 0 and self.cmodel.model_type & 4:
                s = sum(ghmmhelper.arrayint2list(self.cmodel.silent,self.N) )
                if s == 0:
                    # model contains no more silent states
                    self.cmodel.model_type -= 4
            elif silentFlag == 1 and (self.cmodel.model_type & 4) == 0:
                # model contains one silent state
                self.cmodel.model_type += 4


        
        for i in range(self.M):
            ghmmwrapper.set_arrayd(state.b,i,distributionParameters[i])
    
    
    def logprob(self, emissionSequence, stateSequence):
        """ log P[ emissionSequence, stateSequence| m] """
        
        if not isinstance(emissionSequence,EmissionSequence):
            raise TypeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        state = self.getStatePtr( self.cmodel.s, stateSequence[0] )
        emissionProb = ghmmwrapper.get_arrayd(state.b, emissionSequence[0])
        if emissionProb == 0:
            silent = ghmmwrapper.get_arrayint(self.cmodel.silent, 0)
            if silent == 1:
                emissionProb = 1
            else:
                raise SequenceCannotBeBuild, "first symbol " + str(emissionSequence[i+1]) + " not emitted by state " + str(stateSequence[0])
                        
        logP = log(state.pi * emissionProb )
        
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
                            silent = ghmmwrapper.get_arrayint(self.cmodel.silent, stateSequence[i+1] )
                            if silent == 1:
                                emissionProb = 1
                                symbolIndex -= 1
                            else:
                                raise SequenceCannotBeBuild, "symbol " + str(emissionSequence[i+1]) + " not emitted by state "+ str(stateSequence[i+1])

                        logP += log( ghmmwrapper.get_arrayd(state.out_a,j) * emissionProb)
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
            print "Sorry, training of models containing silent states not yet supported."
        else:
            if nrSteps == None:
                ghmmwrapper.reestimate_baum_welch(self.cmodel, trainingSequences.cseq)
            else:
                assert loglikelihoodCutoff != None
                ghmmwrapper.reestimate_baum_welch_nstep(self.cmodel, trainingSequences.cseq,
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
        result = ghmmwrapper.model_apply_background(self.cmodel, cweights)
        
        ghmmwrapper.free_arrayd(cweights)
        if result is not 0:
            print "Doh"
						
    
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
        
        self.cmodel.bp = backgroundObject.cbackground
        self.background = backgroundObject
        self.cmodel.background_id = ghmmhelper.list2arrayint(stateBackground)

        # updating model type
        self.cmodel.model_type += 32 #kHasBackgroundDistributions
    
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
            print str(backgroundID) + " is not contained in labelDomain."    
    
    
    def getBackgroundAssignments(self):
        return ghmmhelper.arrayint2list(self.cmodel.background_id, self.N)
        
    
    def updateTieGroups(self):
        
        assert self.cmodel.tied_to is not None, "cmodel.tied_to is undefined."
        print "##### "+ str(self.cmodel.model_type)
        ghmmwrapper.reestimate_update_tie_groups(self.cmodel)

    
    def setTieGroups(self, tieList):
        
        assert len(tieList) == self.N, "Number of entries in tieList is different from number of states."
        
        if self.cmodel.tied_to is None:
            print "allocating tied_to"
            self.cmodel.tied_to = ghmmhelper.list2arrayint(tieList)
            self.cmodel.model_type += 8
        else:
            print "tied_to already there"
            for i in range(self.N):
                ghmmwrapper.set_arrayint(self.cmodel.tied_to,i,tieList[i])

    def removeTiegroups(self):
        ghmmwrapper.free_arrayi(self.cmodel.tied_to)
        self.cmodel.tied_to = None
        self.model.cmodel.model_type -= 8
    
    def getTieGroups(self):
        assert self.cmodel.tied_to is not None, "cmodel.tied_to is undefined."
        
        return ghmmhelper.arrayint2list(self.cmodel.tied_to, self.N)
    
    def getSilentFlag(self,state):
        return ghmmwrapper.get_arrayint(self.cmodel.silent,state)
    

    def normalize(self):
        """ Normalize transition probs, emission probs (if applicable) """
        
        print "Normalizing now."

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

            print "psum = " + str(pSum)
            for j in range(self.M):
                if pSum >0:  # check for silent state
                    normP = ghmmwrapper.get_arrayd(state.b,j) / pSum
                    ghmmwrapper.set_arrayd(state.b,j,normP)

    def toMatrices(self):
        "Return the parameters in matrix form."
        A = []
        B = []
        pi = []
        for i in range(self.cmodel.N):
            A.append([0.0] * self.N)
            state = self.getStatePtr(self.cmodel.s,i)
            pi.append(state.pi)
            B.append(ghmmhelper.arrayd2list(state.b,self.M ** (state.order+1)))
            for j in range(state.out_states):
                state_index = ghmmwrapper.get_arrayint(state.out_id,j)
                A[i][state_index] = ghmmwrapper.get_arrayd(state.out_a,j)

        return [A,B,pi]


    def toXML(self, filename, backgroundobj = None):
        [A,B,pi] = self.toMatrices()
        nalpha = self.cmodel.M
        nstates = self.cmodel.N

        hmm_dom = xmlutil.HMM()

        hmm_dom.hmmClass.addCode(0, "C1", xmlutil.ValidatingString("Switching class"))

        alphabet = self.emissionDomain
        for c in alphabet.listOfCharacters:
            hmm_dom.hmmAlphabet.addCode(alphabet.index[c], c)

        if backgroundobj != None:
             (n,orders,b) = backgroundobj.tolist()
             background_id = self.getBackgroundAssignments()
             for i in xrange(n):
                  hmm_dom.backgroundDistributions.addDistribution(str(i), orders[i], b[i])

        try:
            tiedlist = self.getTieGroups()
        except:
            print "Ignore tied groups\n"
            print "\"self.cmodel.tiedto\" not defined"
            
        for i in xrange(self.cmodel.N):
            cstate = self.getStatePtr(self.cmodel.s,i)
            if nalpha > 1:
                order = int(log(len(B[i]), nalpha)) - 1
            else:
                order = len(B[i]) - 1

            have_background = False # XXX
            tiedto = None           # XXX
           
            state_dom = xmlutil.HMMState(-1,hmm_dom)
            state_dom.fromDiscreteState( i, pi[i], B[i], cstate.label, order, tiedto, have_background)
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
        assert isinstance(labelDomain,LabelDomain), "Invalid labelDomain"
        self.labelDomain = labelDomain
        
        # Assignment of the C function names to be used with this model type
        self.forwardFunction = ghmmwrapper.foba_logp
        self.forwardAlphaLabelFunction = ghmmwrapper.foba_label_forward
        self.backwardBetaLabelFunction = ghmmwrapper.foba_label_backward
        self.kbestFunction = ghmmwrapper.kbest        
        self.gradientDescentFunction = ghmmwrapper.gradient_descent
        self.cmodel.state_label = ghmmwrapper.int_array(self.N)

    def __str__(self):
        hmm = self.cmodel
        strout = "\nGHMM Model\n"
        strout += "Name: " + str(self.cmodel.name)
        strout += "\nModelflags: "+ self.printtypes(self.cmodel.model_type)
        strout += "\nNumber of states: "+ str(hmm.N)
        strout += "\nSize of Alphabet: "+ str(hmm.M)
        for k in range(hmm.N):
            state = ghmmwrapper.get_stateptr(hmm.s, k)
            strout += "\n\nState number "+ str(k) +":"

            #strout += "\nState label: "+ str(self.r_index[state.label])
            strout += "\nState label: "+ str(self.labelDomain.external(state.label))

            strout += "\nState order: " + str(state.order)
            strout += "\nInitial probability: " + str(state.pi)
            #strout += "\nsilent state: " + str(get_arrayint(self.model.silent,k))
            strout += "\nOutput probabilites:\n"
            for outp in range(hmm.M**(state.order+1)):
                strout+=str(ghmmwrapper.get_arrayd(state.b,outp))
                if outp % hmm.M == hmm.M-1:
                    strout += "\n"
                else:
                    strout += ", "

            strout += "Outgoing transitions:"
            for i in range( state.out_states):
                strout += "\ntransition to state " + str( ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.out_a,i))
            strout +=  "\nIngoing transitions:"
            for i in range(state.in_states):
                strout +=  "\ntransition from state " + str( ghmmwrapper.get_arrayint(state.in_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.in_a,i))
                strout += "\nint fix:" + str(state.fix) + "\n"
        strout += "\nSilent states: \n"
        for k in range(hmm.N):
            strout += str(ghmmwrapper.get_arrayint(self.cmodel.silent,k)) + ", "
        strout += "\n"
        return strout

    def setLabels(self,labelList):
        """  Set the state labels to the values given in labelList.
             LabelList is in external representation.
        """
        
        assert len(labelList) == self.N, "Invalid number of labels."
        
        # set state label to to the appropiate index
        for i in range(self.N):
            if not self.labelDomain.isAdmissable(labelList[i]):
                raise GHMMOutOfDomain, "Label "+str(labelList[i])+" not included in labelDomain."
            state = self.getStatePtr(self.cmodel.s, i)
            state.label = self.labelDomain.internal(labelList[i])

    def getLabels(self):
        labels = []
        for i in range(self.N):
            state = self.getStatePtr(self.cmodel.s, i)
            labels.append(self.labelDomain.external(state.label))
        
        return labels    
    
    def getLabel(self,stateIndex):
        """ Returns label of the state 'stateIndex'.
        
        """
        state = self.getStatePtr(self.cmodel.s, stateIndex)
        return state.label
         
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
        seqPtr = ghmmwrapper.model_label_generate_sequences(self.cmodel,seed,seqLength,1,seqLength)
        return EmissionSequence(self.emissionDomain, seqPtr, labelDomain = self.labelDomain )

    def sample(self, seqNr,seqLength, seed = 0):
        seqPtr = ghmmwrapper.model_label_generate_sequences(self.cmodel,seed,seqLength,seqNr,seqLength)
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
            print "Sorry, k-best decoding on models containing silent states not yet supported."
        else:
            if isinstance(emissionSequences,EmissionSequence):
                seqNumber = 1
            elif isinstance(emissionSequences,SequenceSet):
                seqNumber = len(emissionSequences)
            else:
                raise TypeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)

            log_p = ghmmwrapper.double_array(0)

            allLogs = []
            allLabels = []
            seq_todo = seqNumber
            print " "
            for i in range(seqNumber):
                if seqNumber > 1:
                    print seq_todo
                    seq_todo -= 1
                seq = emissionSequences.getPtr(emissionSequences.cseq.seq, i)
                seq_len = ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i)

                labeling =  self.kbestFunction(self.cmodel, seq, seq_len, k, log_p)
                oneLabel = []

                for j in range(seq_len):
                    oneLabel.append(ghmmwrapper.get_arrayint(labeling,j))

                allLabels.append(oneLabel)
                allLogs.append(ghmmwrapper.get_arrayd(log_p, 0))
                ghmmwrapper.free_arrayi(labeling)

                            
            ghmmwrapper.free_arrayd(log_p)
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


        cmodelPTR = ghmmwrapper.discrime_modelarray_alloc(1)
        ghmmwrapper.discrime_modelarray_setptr(cmodelPTR, self.cmodel, 0)
        error = self.gradientDescentFunction(cmodelPTR, emissionsequences.cseq, eta, steps)

        self.cmodel = ghmmwrapper.discrime_modelarray_getptr(cmodelPTR, 0)
        ghmmwrapper.discrime_modelarray_dealloc(cmodelPTR)
        
        if error == -1:
            print "ERROR: Gradient descent finished not successfully."

        return error
        
    def modelNormalize(self):
       i_error = ghmmwrapper.model_normalize(self.cmodel)
       if i_error == -1:
            print "ERROR: normalize finished with -1"
       

        
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
            ret_val = ghmmwrapper.foba_label_logp(self.cmodel, seq, labels, tmp, likelihood)

            if ret_val == -1:

                print "Warning: forward returned -1: Sequence", i,"cannot be build."
                likelihoodList.append(-float('Inf'))
            else:
                likelihoodList.append(ghmmwrapper.get_arrayd(likelihood,0))

        ghmmwrapper.free_arrayd(likelihood)
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

        calpha = ghmmwrapper.matrix_d_alloc(t,n_states)
        cscale = ghmmwrapper.double_array(t)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)
        label = ghmmwrapper.int_array(t)

        for i in range(len(labelSequence)):
            ghmmwrapper.set_arrayint(label, i, self.internalLabel(labelSequence[i]))

        error = self.forwardAlphaLabelFunction(self.cmodel, seq, label, t, calpha, cscale, logP)
        if error == -1:
            print "ERROR: Forward finished with -1: Sequence " + str(i) + " cannot be build."
           

        # translate alpha / scale to python lists
        pyscale = ghmmhelper.arrayd2list(cscale, t)
        pyalpha = ghmmhelper.matrixd2list(calpha,t,n_states)
        logpval = ghmmwrapper.get_arrayd(logP, 0)
        
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
            print "ERROR: backward finished with -1: EmissionSequence cannot be build."
            

        pybeta = ghmmhelper.matrixd2list(cbeta,t,self.cmodel.N)
        logpval = ghmmwrapper.get_arrayd(logP, 0)
        
        # deallocation
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
            print "Sorry, training of models containing silent states not yet supported."
        else:
            if nrSteps == None:
                ghmmwrapper.reestimate_baum_welch_label(self.cmodel, trainingSequences.cseq)
            else:
                assert loglikelihoodCutoff != None
                ghmmwrapper.reestimate_baum_welch_nstep_label(self.cmodel, trainingSequences.cseq,
                                                        nrSteps, loglikelihoodCutoff)


class GaussianEmissionHMM(HMM):
    """ HMMs with Gaussian distribution as emissions.
        
    """

    def __init__(self, emissionDomain, distribution, cmodel):
        HMM.__init__(self, emissionDomain, distribution, cmodel)
        self.silent = 0  # flag indicating whether the model type does include silent states

        # Assignment of the C function names to be used with this model type
        self.freeFunction = ghmmwrapper.call_smodel_free
        self.samplingFunction = ghmmwrapper.smodel_generate_sequences
        self.viterbiFunction = ghmmwrapper.sviterbi
        self.forwardFunction = ghmmwrapper.sfoba_logp
        self.forwardAlphaFunction = ghmmwrapper.sfoba_forward
        self.backwardBetaFunction = ghmmwrapper.sfoba_backward
        self.getStatePtr = ghmmwrapper.get_sstate_ptr
        self.fileWriteFunction = ghmmwrapper.call_smodel_print
        self.getModelPtr = ghmmwrapper.get_smodel_ptr
        self.castModelPtr = ghmmwrapper.cast_smodel_ptr
        self.distanceFunction = ghmmwrapper.smodel_prob_distance

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
    
    def __str__(self):
        hmm = self.cmodel
        strout = "\nHMM Overview:"
        strout += "\nNumber of states: " + str(hmm.N)
        strout += "\nNumber of mixture components: " + str(hmm.M)

        for k in range(hmm.N):
            state = ghmmwrapper.get_sstate(hmm, k)
            strout += "\n\nState number "+ str(k) + ":"
            strout += "\nInitial probability: " + str(state.pi) + "\n"

            weight = ""
            mue = ""
            u =  ""

            weight += str(ghmmwrapper.get_arrayd(state.c,0))
            mue += str(ghmmwrapper.get_arrayd(state.mue,0))
            u += str(ghmmwrapper.get_arrayd(state.u,0))

            strout += "  mean: " + str(mue) + "\n"
            strout += "  variance: " + str(u) + "\n"
            strout += "  fix: " + str(state.fix) + "\n"
            
            for c in range(self.cmodel.cos):
                strout += "\n  Class : " + str(c)                
                strout += "\n    Outgoing transitions:"
                for i in range( state.out_states):
                    strout += "\n      transition to state " + str(ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability = " + str(ghmmwrapper.get_2d_arrayd(state.out_a,c,i))
                strout +=  "\n    Ingoing transitions:"
                for i in range(state.in_states):
                    strout += "\n      transition from state " + str(ghmmwrapper.get_arrayint(state.in_id,i) ) +" with probability = "+ str(ghmmwrapper.get_2d_arrayd(state.in_a,c,i))


        return strout

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
        calpha = ghmmwrapper.matrix_d_alloc(t,i)
        cscale = ghmmwrapper.double_array(t)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)

        error = self.forwardAlphaFunction(self.cmodel, seq,t, None, calpha, cscale, logP)
        if error == -1:
            print "ERROR: Forward finished with -1: Sequence " + str(seq_nr) + " cannot be build."
            

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
            print "ERROR: backward finished with -1: EmissionSequence cannot be build."
            

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
            print "self.cmodel.cos = ", self.cmodel.cos
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
                
                print "Warning: forward returned -1: Sequence", i,"cannot be build."
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

        ghmmwrapper.free_arrayd(likelihood)
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
            print "self.cmodel.cos = ", self.cmodel.cos
            assert self.cmodel.class_change is not None, "Error: class_change not initialized."


        log_p = ghmmwrapper.double_array(1)

        allLogs = []
        allPaths = []
        for i in range(seqNumber):
            if self.cmodel.cos > 1:
                self.cmodel.class_change.k = i

            seq = emissionSequences.getPtr(emissionSequences.cseq.seq,i)
            seq_len = ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i)

            if seq_len != 0:
                viterbiPath = self.viterbiFunction(self.cmodel,seq,seq_len,log_p)
            else:
                viterbiPath = None

            onePath = []
            
            # for model types without possible silent states the length of the viterbi path is known
            if self.silent == 0:            
                for j in range(seq_len):                
                    onePath.append(ghmmwrapper.get_arrayint(viterbiPath,j))
            
            # in the silent case we have to reversely append as long as the path is positive because unused positions
            # are initialised with -1 on the C level.
            elif self.silent == 1:
                
                for j in range( ( seq_len * self.N )-1,-1,-1): # maximum length of a viterbi path for a silent model
                    d = ghmmwrapper.get_arrayint(viterbiPath,j)
                                   
                    if d >= 0:
                        onePath.insert(0,d)
                    else:
                        break

            allPaths.append(onePath)
            allLogs.append(ghmmwrapper.get_arrayd(log_p, 0))
        
        ghmmwrapper.free_arrayd(log_p)
        ghmmwrapper.free_arrayi(viterbiPath)  
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
            raise TypeError, "baumWelch requires a SequenceSet object."
        
        assert self.emissionDomain.CDataType == "double", "Continous sequence needed."
        
        self.baumWelchSetup(trainingSequences, nrSteps)
        ghmmwrapper.sreestimate_baum_welch(self.BWcontext)        
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

            Return: a C-array of type smosqd_t
        """
        self.BWcontext = ghmmwrapper.smosqd_t_array(1)
        self.BWcontext.smo  = self.cmodel
        self.BWcontext.sqd  = trainingSequences.cseq    # copy reference to sequence_d_t
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

        ghmmwrapper.free_arrayd(self.BWcontext.logp)
        self.BWcontext.logp = None
        ghmmwrapper.free_smosqd_t(self.BWcontext)
        self.BWcontext = None


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

    def __str__(self):
        "defines string representation"
        hmm = self.cmodel

        strout = "\nOverview of HMM:"
        strout += "\nNumber of states: "+ str(hmm.N)
        strout += "\nNumber of output distributions per state: "+ str(hmm.M)

        for k in range(hmm.N):
            state = ghmmwrapper.get_sstate(hmm, k)
            strout += "\n\nState number "+ str(k) +":"
            strout += "\nInitial probability: " + str(state.pi)
            strout += "\n"+ str(hmm.M) + " density function(s):\n"

            weight = ""
            mue = ""
            u =  ""

            for outp in range(hmm.M):
                weight += str(ghmmwrapper.get_arrayd(state.c,outp))+", "
                mue += str(ghmmwrapper.get_arrayd(state.mue,outp))+", "
                u += str(ghmmwrapper.get_arrayd(state.u,outp))+", "

            strout += "  pdf component weights : " + str(weight) + "\n"
            strout += "  mean vector: " + str(mue) + "\n"
            strout += "  variance vector: " + str(u) + "\n"
            strout += "\nOutgoing transitions:"

            for i in range( state.out_states):
                strout += "\ntransition to node " + str( ghmmwrapper.get_arrayint(state.out_id,i) ) +" with probablity = "+ str(ghmmwrapper.get_2d_arrayd(state.out_a,0,i))


            strout +=  "\nIngoing transitions:"

            for i in range(state.in_states):
                strout += "\ntransition from node " + str( ghmmwrapper.get_arrayint(state.in_id,i) ) + " with probablity = "+ str(ghmmwrapper.get_2d_arrayd(state.in_a,0,i))

            strout += "\nint fix:" + str(state.fix) + "\n"
        return strout

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


def HMMDiscriminativeTraining(HMMList, SeqList, gradient = 0):
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

    HMMArray = ghmmwrapper.discrime_modelarray_alloc(inplen)
    SeqArray = ghmmwrapper.discrime_seqarray_alloc(inplen)

    for i in range(inplen):
        ghmmwrapper.discrime_modelarray_setptr(HMMArray, HMMList[i].cmodel, i)
        ghmmwrapper.discrime_seqarray_setptr(SeqArray, SeqList[i].cseq, i)

    ghmmwrapper.discriminative(HMMArray, SeqArray, inplen, gradient)

    for i in range(inplen):
        HMMList[i].cmodel = ghmmwrapper.discrime_modelarray_getptr(HMMArray, i)
        SeqList[i].cseq   = ghmmwrapper.discrime_seqarray_getptr(SeqArray, i)

    ghmmwrapper.discrime_modelarray_dealloc(HMMArray)
    ghmmwrapper.discrime_seqarray_dealloc(SeqArray)

    return HMMDiscriminativePerformance(HMMList, SeqList)


        
        
def HMMDiscriminativePerformance(HMMList, SeqList):

    if len(HMMList) != len(SeqList):
        raise TypeRrror, 'Inputs not equally long'

    inplen = len(HMMList)
    
    single = [0.0] * inplen

    HMMArray = ghmmwrapper.discrime_modelarray_alloc(inplen)
    SeqArray = ghmmwrapper.discrime_seqarray_alloc(inplen)

    for i in range(inplen):
        ghmmwrapper.discrime_modelarray_setptr(HMMArray, HMMList[i].cmodel, i)
        ghmmwrapper.discrime_seqarray_setptr(SeqArray, SeqList[i].cseq, i)

    retval = ghmmwrapper.discrime_compute_performance(HMMArray, SeqArray, inplen)
    
    ghmmwrapper.discrime_modelarray_dealloc(HMMArray)
    ghmmwrapper.discrime_seqarray_dealloc(SeqArray)

    return retval
        
       
        
