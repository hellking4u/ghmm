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

   # XXX do we need distribution here? I think we need it to interpret the parameters in B.

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
import re
from os import path

ghmmwrapper.gsl_rng_init() #Init random number generator

#-------------------------------------------------------------------------------
#- Exceptions ------------------------------------------------------------------

class GHMMError(Exception):
    """Base class for exceptions in this module."""
    def __init__(self, message):
        self.message = message

class UnknownInputType(GHMMError):
    def __init__(self,message):
       print "\n\n UnknownInputType Exception: " + str(message) + "\n"

class NoValidCDataType(GHMMError):
    def __init__(self,message):   
        print "\n\n NoValidCDataType: " + str(message) + "\n"
        
class badCPointer(GHMMError):
    def __init__(self,message):   
        print "\n\nbadCPointer Exception: " + str(message) + "\n"
                
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


DNA = Alphabet(['a','c','g','t'])
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
        pass
        #return isinstance(emission,Float)
   

    
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
        self.prob_vector = None
        
    def set(self, prob_vector):
        self.prob_vector = prob_vector

    def get(self):
        return self.prob_vector
        

class ContinousDistribution(Distribution):
    pass

class GaussianDistribution(ContinousDistribution):

    def __init__(self, domain):
        self.emissionDomain = domain
        self.mu = None
        self.sigma = None
        
    def set(self, (mu, sigma)):
        self.mu = mu
        self.sigma = sigma

    def get(self):
        return (self.mu, self.sigma2)
        

class MixtureContinousDistribution(Distribution):
    pass


#-------------------------------------------------------------------------------
#- Sequence, SequenceSet and derived  ------------------------------------------
class EmissionSequence:
    """ An EmissionSequence contains the *internal* representation of
        a sequence of emissions. It also contains a reference to the
        domain where the emission orginated from.
    """

    def __init__(self, emissionDomain, sequenceInput):
        self.emissionDomain = emissionDomain 
        # XXX How do you maintain reference to the SequenceSet.cseq?
        
        if self.emissionDomain.CDataType == "int": # underlying C data type is integer
            self.getPtr = ghmmwrapper.get_col_pointer_int # defines C function to be used to access a single sequence
            self.getSymbol = ghmmwrapper.get_2d_arrayint
            self.setSymbol = ghmmwrapper.set_2d_arrayint
            
            if isinstance(sequenceInput, list):  
                (seq,l) = ghmmhelper.list2matrixint([sequenceInput])
                self.cseq = ghmmwrapper.sequence_calloc(1)
                self.cseq.seq = seq
                self.cseq.seq_number = 1
                ghmmwrapper.set_arrayint(self.cseq.seq_len,0,l[0]) 
                
            elif isinstance(sequenceInput, str): # from file
                # Does it make sense to read only one sequence from file?? 
                # sequence_t **sequence_read(const char *filename, int *seq_arrays)
                pass
            elif isinstance(sequenceInput, ghmmwrapper.sequence_t):# internal use
                # XXX Should maintain reference to the parent
                if sequenceInput.seq_number > 1:
                    raise badCPointer, "Use SequenceSet for multiple sequences."
                self.cseq = sequenceInput
                
            else:    
                raise UnknownInputType, "inputType " + str(type(sequenceInput)) + " not recognized."
        
        elif self.emissionDomain.CDataType == "double": # underlying C data type is double
            self.getPtr = ghmmwrapper.get_col_pointer_d # defines C function to be used to access a single sequence
            self.getSymbol = ghmmwrapper.get_2d_arrayd
            self.setSymbol = ghmmwrapper.set_2d_arrayd
                        
            if isinstance(sequenceInput, list): 
                (seq,l) = ghmmhelper. list2matrixd([sequenceInput])
                self.cseq = ghmmwrapper.sequence_d_calloc(1)
                self.cseq.seq = seq
                self.cseq.seq_number = 1
                set_arrayint(self.cseq.seq_len,0,1)
				

            elif isinstance(sequenceInput, str): # from file
                pass
            elif isinstance(sequenceInput, ghmmwrapper.sequence_d_t): # internal use
                # XXX Should maintain reference to the parent
                if sequenceInput.seq_number > 1:
                    raise badCPointer, "Use SequenceSet for multiple sequences."
                self.cseq = sequenceInput
            else:    
                raise UnknownInputType, "inputType " + str(type(sequenceInput)) + " not recognized."
        
        
        else:
            raise NoValidCDataType, "C data type " + str(self.emissionDomain.CDataType) + " invalid."

        
    def __del__(self):        
        if self.emissionDomain.CDataType == "int": # delete the C ptr to sequence_t
            ghmmwrapper.sequence_clean(self.cseq)
            
        elif self.emissionDomain.CDataType == "double": # delete the C ptr to sequence_d_t
            ghmmwrapper.sequence_d_clean(self.cseq)

        
    def __len__(self):
        "In order to be consistent with SequenceSet len(self) returns the number of sequences."
        return ghmmwrapper.get_arrayint(self.cseq.seq_len,0)	


    def __setitem__(self, index, value):
        self.setSymbol(self.seq_c.seq,0,index,value)
		

    def __getitem__(self, index):
        """ Return the symbol at position 'index'. """
        if index < len(self):
            return self.getSymbol(self.cseq.seq, 0, index)
        
        else:
            raise IndexError    
    
    def __str__(self):
        "Defines string representation."
        seq = self.cseq
        strout = ""
        strout += "\nEmissionSequence Instance:\nlength " + str(ghmmwrapper.get_arrayint(seq.seq_len,0))+ ", weight " + str(ghmmwrapper.get_arrayd(seq.seq_w,0))  + ":\n"
        for j in range(ghmmwrapper.get_arrayint(seq.seq_len,0) ):
            strout += str( self.emissionDomain.external(self[j]) )   
           
    	return strout		

#    def sequenceSet(self):
#        """ Make a one-element SequenceSet out of me """


class SequenceSet:
    def __init__(self, emissionDomain, sequenceSetInput):
        self.emissionDomain = emissionDomain 
        self.cseq = None
        
        if self.emissionDomain.CDataType == "int": # underlying C data type is integer
            # necessary C functions for accessing the sequence_t struct
            self.calloc = ghmmwrapper.sequence_calloc
            self.getPtr = ghmmwrapper.get_col_pointer_int # defines C function to be used to access a single sequence
            self.castPtr = ghmmwrapper.cast_ptr_int # cast int* to int** pointer
            self.emptySeq = ghmmwrapper.int_2d_array_nocols # allocate an empty int**
            self.setSeq = ghmmwrapper.set_2d_arrayint_col # assign an int* to a position within a int**
            
            if isinstance(sequenceSetInput, str): # from file
                # reads in the first sequence struct in the input file
                self.cseq  = ghmmwrapper.seq_read(sequenceSetInput)
                
            elif isinstance(sequenceSetInput, list): 
                seq_nr = len(sequenceSetInput)
                self.cseq = ghmmwrapper.sequence_calloc(seq_nr)
                self.cseq.seq_number = seq_nr

                (seq,lenghts) = ghmmhelper.list2matrixint(sequenceSetInput)

                self.cseq.seq = seq
                for i in range(seq_nr):
                    ghmmwrapper.set_arrayint(self.cseq.seq_len,i ,lenghts[i])


            elif isinstance(sequenceSetInput, ghmmwrapper.sequence_t): # inputType == sequence_t**
                self.cseq = sequenceSetInput
                
            else:    
                raise UnknownInputType, "inputType " + str(type(sequenceSetInput)) + " not recognized."
    
            self.__array = []    
            for index in range(len(self)):
                oneseq = self.getPtr(self.cseq.seq,index)
                self.__array.append(oneseq)    
        
        elif self.emissionDomain.CDataType == "double": # underlying C data type is double
            # necessary C functions for accessing the sequence_d_t struct
            self.calloc =  ghmmwrapper.sequence_d_calloc
            self.getPtr = ghmmwrapper.get_col_pointer_d # defines C function to be used to access a single sequence
            self.castPtr = ghmmwrapper.cast_ptr_d # cast double* to double** pointer
            self.emptySeq = ghmmwrapper.double_2d_array_nocols  # cast double* to int** pointer
            self.setSeq = ghmmwrapper.set_2d_arrayd_col # assign an double* to a position within a double**
                        
            if isinstance(sequenceSetInput, list): 
                seq_nr = len(sequenceSetInput)
                self.cseq = ghmmwrapper.sequence_d_calloc(seq_nr)
                self.cseq.seq_number = seq_nr

                (seq,lenghts) = ghmmhelper.list2matrixd(sequenceSetInput)
                self.cseq.seq = seq
                for i in range(seq_nr):
                    ghmmwrapper.set_arrayint(self.cseq.seq_len, i, lenghts[i])

            elif isinstance(sequenceSetInput, str): # from file
                print "fromFile", sequenceSetInput                
                self.cseq  = ghmmwrapper.seq_d_read(sequenceSetInput)
                                
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
        if self.emissionDomain.CDataType == "int":
            ghmmwrapper.sequence_clean(self.cseq)
            
        elif self.emissionDomain.CDataType == "double":
            ghmmwrapper.sequence_d_clean(self.cseq)
            
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
        

    def weight(self, i):
        """ Return the weight of sequence i. Weights are used in Baum-Welch"""
        return ghmmwrapper.get_arrayd(self.cseq.seq_w,i)

    def setWeight(self, i, w):
        """ Return the weight of sequence i. Weights are used in Baum-Welch"""
        return ghmmwrapper.set_arrayd(self.cseq.seq_w,i,w)
        
    
    def __setitem__(self, index, value): # Only allow EmissionSequence?
        pass
        
    def __getitem__(self, index):
        """ Return an EmissionSequence object initialized with sequence 'index'.
        
        """
        
        # Right now only the pointer is passed. Do we want to have a real copy of the sequence instead ?
        
        seq = self.calloc(1)
        seq.seq = self.castPtr(self.__array[index]) # int* -> int** reference
        ghmmwrapper.set_arrayint(seq.seq_len,0,ghmmwrapper.get_arrayint(self.cseq.seq_len,index))
        seq.seq_number = 1
        
        return EmissionSequence(self.emissionDomain, seq)

    def getSubset(self, seqIndixes):
        """ Returns a SequenceSet containing (references to) the sequences with the indixes in
            'seqIndixes'.
       
        """
        seqNumber = len(seqIndixes)       
        seq = ghmmwrapper.sequence_calloc(seqNumber)
        seq.seq = self.emptySeq(seqNumber)
        for i in range(seqNumber):
            self.setSeq(seq.seq,i,self.__array[ seqIndixes[i] ]) 
            ghmmwrapper.set_arrayint(seq.seq_len,i,ghmmwrapper.get_arrayint(self.cseq.seq_len,seqIndixes[i]))
        seq.seq_number = seqNumber
        
        return SequenceSet(self.emissionDomain, seq)

        
def SequenceSetOpen(emissionDomain, fileName):
    """ Reads a sequence file with multiple sequence sets. 
    
    Returns a list of SequenceSet objects.
    
    """
    
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
       
    ghmmwrapper.freearray(dArr)    
    return  sequenceSets



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

    def __call__(self, fileName, modelIndex = None):
        # MO & SMO Files
        (hmmClass, emission_domain, distribution) = self.determineHMMClass(fileName)
        nrModelPtr = ghmmwrapper.int_array(1)
        models = ghmmwrapper.smodel_read(fileName, nrModelPtr)
        nrModels = ghmmwrapper.get_arrayint(nrModelPtr, 0)
        if modelIndex == None:
            cmodel = get_smodel_ptr(models, 0) # XXX Who owns the pointer?
        else:
            if modelIndex < nrModels:
                cmodel = get_smodel_ptr(models, modelIndex) # XXX Who owns the pointer?
            else:
                print "modelIndex too large -- only have ", nrModels
                return None
        m = hmmClass(emission_domain, distribution(emission_domain), cmodel)
        return m


    def all(self, fileName):
        # MO & SMO Files
        (hmmClass, emission_domain, distribution) = self.determineHMMClass(fileName)
        nrModelPtr = ghmmwrapper.int_array(1)
        models = ghmmwrapper.smodel_read(fileName, nrModelPtr)
        nrModels = ghmmwrapper.get_arrayint(nrModelPtr, 0)
        result = []
        for i in range(nrModels):
            cmodel = ghmmwrapper.get_smodel_ptr(models, i)
            m = hmmClass(emission_domain, distribution(emission_domain), cmodel)
            result.append(m)
        return result
        

    def determineHMMClass(self, fileName):
        #
        # smo files 
        #
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

        return (None, None, None)
            

        

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

                silent_flag = 0
                silent_states = []

                #initialize states
                for i in range(cmodel.N):
                    state = ghmmwrapper.get_stateptr(states,i)
                    state.b = ghmmhelper.list2arrayd(B[i])
                    state.pi = pi[i]
                    
                    if (sum(B[i]) == 0 ): 
                        silent_states.append(1)
                        silent_flag = 4
                    else:
                        silent_states.append(0)

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
                cmodel.model_type = silent_flag
                cmodel.silent = ghmmhelper.list2arrayint(silent_states)
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

                    print "  state " + str(i) + ":"
                    
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
        
        self.N = self.cmodel.N  # number of states
        self.M = self.cmodel.M  # number of symbols / mixture components
        
        self.silent = 0   # flag for silent states
        
        # Function pointer to the C level (defined in derived classes)
        self.freeFunction = ""           # C function for deallocation
        self.samplingFunction = ""       # C function for sequence generation
        self.viterbiFunction = ""        # C function viterbi algorithm
        self.forwardFunction = ""        # C function forward algortihm (likelihood only)
        self.forwardAlphaFunction = ""   # C function forward algorithm (alpha matrix,scale)
        self.backwardBetaFunction = ""   # C function backkward algorithm (beta matrix)
        self.getStatePtr = ""            # C function to get a pointer to a state struct
       # self.getModelPtr = ""            # C function to get a pointer to the model struct
        
        
  #  def __del__(self):
  #      """ Deallocation routine for the underlying C data structures. """
  #      modelPtr = getModelPtr(self.cmodel)
  #      self.freeFunction(modelPtr)
        
    
    
    def loglikelihood(self, emissionSequences): 
        """ Compute log( P[emissionSequences| model]) using the forward algorithm
            assuming independence of the sequences in emissionSequences

            emissionSequences can either be a SequenceSet or a Sequence

            Result: log( P[emissionSequences| model]) of type float which is
            computed as \sum_{s} log( P[s| model]) when emissionSequences
            is a SequenceSet
            
            Note: The implementation will not compute the full forward matrix (XXX ToDo)
        """
        
        if not isinstance(emissionSequences,EmissionSequence) and not isinstance(emissionSequences,SequenceSet):
            raise AttributeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)        
        

        
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
            raise AttributeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)        
              
        likelihood = ghmmwrapper.double_array(1)
        likelihoodList = []
        
        
        for i in range(seqNumber):
            seq = emissionSequences.getPtr(emissionSequences.cseq.seq,i)
            ret_val = self.forwardFunction(self.cmodel, seq,ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i), likelihood)
            likelihoodList.append(ghmmwrapper.get_arrayd(likelihood,0))
            # XXX check ret_val

        ghmmwrapper.free_arrayd(likelihood)  
        return likelihoodList

       

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

        if not isinstance(trainingSequences, SequenceSet):
            raise TypeError
            
        if (self.emissionDomain.CDataType == "double"):
            smosqd_cpt = self.baumWelchSetup(trainingSequences, nrSteps)
            ghmmwrapper.sreestimate_baum_welch(smosqd_cpt)        
            likelihood = ghmmwrapper.get_arrayd(smosqd_cpt.logp, 0)
            ghmmwrapper.free_smosqd_t(smosqd_cpt)
            
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
        baumWelchCData = ghmmwrapper.smosqd_t_array(1)

        smosqd_ptr = ghmmwrapper.get_smosqd_t_ptr(baumWelchCData, 0)
        smosqd_ptr.smo  = self.cmodel
        smosqd_ptr.sqd  = trainingSequences.cseq    # copy reference to sequence_d_t
        smosqd_ptr.logp = ghmmwrapper.double_array(1) # place holder for sum of log-likelihood
        smosqd_ptr.eps  = 10e-6
        smosqd_ptr.max_iter = nrSteps
           
        return baumWelchCData
    
    def baumWelchStep(self, nrSteps, loglikelihoodCutoff):
        """ Compute one iteration of Baum Welch estimation.
            Use baum_welch_setup and baum_welch_step if you want more control
            over the training, compute diagnostics or do noise-insertion

            training_sequences can either be a SequenceSet or a Sequence
        """
        pass
    
    def baumWelchDelete(self):
        """ Delete the necessary temporary variables for Baum-Welch-reestimation """
        pass
        
    def forward(self, emissionSequence):
        """

            Result: the (N x T)-matrix containing the forward-variables
                    and the scaling vector
        """
                      
        if not isinstance(emissionSequence,EmissionSequence):
            raise AttributeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        unused = ghmmwrapper.double_array(1) # Dummy return value for forwardAlphaFunction    	
     
        t = len(emissionSequence)
        calpha = ghmmwrapper.matrix_d_alloc(t,self.N)
        cscale = ghmmwrapper.double_array(t)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)	
		
        error = self.forwardAlphaFunction(self.cmodel, seq,t, calpha, cscale, unused)
        if error == -1:
            print "ERROR: forward finished with -1: EmissionSequence cannot be build."
            exit
        
        # translate alpha / scale to python lists 
        pyscale = ghmmhelper.arrayd2list(cscale, t)
        pyalpha = ghmmhelper.matrixd2list(calpha,t,self.N)
        
        # deallocation
        ghmmwrapper.freearray(unused)
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(calpha,t)
        return (pyalpha,pyscale) 
        
        

    def backward(self, emissionSequence, scalingVector):
        """

            Result: the (N x T)-matrix containing the backward-variables
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise AttributeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)
        
        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2arrayd(scalingVector)
        
        # alllocating beta matrix
        t = len(emissionSequence)
        cbeta = ghmmwrapper.double_2d_array(t, self.N)
        
        error = self.backwardBetaFunction(self.cmodel,seq,t,cbeta,cscale)
        if error == -1:
            print "ERROR: backward finished with -1: EmissionSequence cannot be build."
            exit
        
        pybeta = ghmmhelper.matrixd2list(cbeta,t,self.N)

        # deallocation
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(cbeta,t)
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
            raise AttributeError, "EmissionSequence or SequenceSet required, got " + str(emissionSequences.__class__.__name__)                            


        log_p = ghmmwrapper.double_array(1)
        
        allPaths = []
        for i in range(seqNumber):
            seq = emissionSequences.getPtr(emissionSequences.cseq.seq,i)
            seq_len = ghmmwrapper.get_arrayint(emissionSequences.cseq.seq_len,i)
            viterbiPath =  self.viterbiFunction(self.cmodel,seq,seq_len,log_p)

            #viterbi_prob = get_arrayd( log_p, 0 )
        
            onePath = []
            
            # for model types without possible silent states the length of the viterbi path is known
            if self.silent == 0:            
                for i in range(seq_len):                
                    onePath.append(ghmmwrapper.get_arrayint(viterbiPath,i))
            
            # in the silent case we have to reversely append as long as the path is positive because unused positions
            # are initialised with -1 on the C level.
            elif self.silent == 1:   
                
                for i in range( ( seq_len * self.N )-1,-1,-1): # maximum length of a viterbi path for a silent model
                    d = ghmmwrapper.get_arrayint(viterbiPath,i)
                                   
                    if d >= 0:
                        onePath.insert(0,d)
                    else:
                        break
                                        
            allPaths.append(onePath)
        if emissionSequences.cseq.seq_number > 1:
            return allPaths
        else:
            return allPaths[0]    
                    
    
    def sample(self, seqNr ,T):
        """ Sample emission sequences 


        """
        seqPtr = self.samplingFunction(self.cmodel,0,T,seqNr,self.N)
        return SequenceSet(self.emissionDomain,seqPtr)
        

    def sampleSingle(self, T):
        """ Sample a single emission sequence of length at most T.
            Returns a Sequence object.
        """
        seqPtr = self.samplingFunction(self.cmodel,0,T,1,self.N)
        return EmissionSequence(self.emissionDomain,seqPtr)

    def state(self, stateLabel):
        """ Given a stateLabel return the integer index to the state """
        pass

    def getInitial(self, i):
        """ Accessor function for the initial probability \pi_i """
        return ghmmwrapper.get_arrayd(self.cmodel.pi,i)

    def setInitial(self, i, prob, fixProb=0):
        """ Accessor function for the initial probability \pi_i
            For 'fixProb' = 1 \pi will be rescaled to 1 with 'prob' fixed to the
            arguement value.
           
         """
        ghmmwrapper.set_array_d(self.cmodel.pi,i,j)
        # renormalizing pi, pi(i) is fixed on value 'prob'
        if fixProb == 1:
            coeff = 1 - i
            for j in range(self.N):
                if i != j:
                    pi_j = ghmmwrapper.get_arrayd(self.cmodel.pi,j)    
                    pi_j = pi_j * i / coeff
                    ghmmwrapper.set_array_d(self.cmodel.pi,j,pi_j)

    def getTransition(self, i, j):
        """ Accessor function for the transition a_ij """
        state = self.getStatePtr(self.cmodel.s,i)
        transition = 0
        for i in range(state.out_states):
            stateId = ghmmwrapper.get_arrayint(state.out_id,i)
            if stateId == j:
                transition = ghmmwrapper.get_arrayd(state.out_a,i)
                break
        return transition    
            
    def setTransition(self, i, j, prob):
        """ Accessor function for the transition a_ij """
        pass

    def getEmission(self, i):
        """ Accessor function for the emission distribution parameters of state 'i'. 
        
            For discrete models the distribution over the symbols is returned,
            for continous models a matrix of the form 
            [ [mu_1, sigma_1, weight_1] ... [mu_M, sigma_M, weight_M]  ] is returned.
                    
        """
        if self.emissionDomain.CDataType == "int": # discrete emissions.
            state = self.getStatePtr(self.cmodel.s,i)
            emissions = ghmmhelper.arrayd2list(state.b,self.M)
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


def HMMwriteList(fileName,hmmList):
    if path.exists(fileName):
        print "HMMwriteList warning: File " + str(fileName) + " already exists. " + str(len(hmmList)) + " new models will be appended."
    for model in hmmList:
        model.write(fileName)
        
    

        

class DiscreteEmissionHMM(HMM):
    
    def __init__(self, emissionDomain, distribution, cmodel):
        HMM.__init__(self, emissionDomain, distribution, cmodel)
        self.silent = 1  # flag indicating whether the model type does include silent states
        self.samplingFunction = ghmmwrapper.model_generate_sequences
        self.viterbiFunction = ghmmwrapper.viterbi
        self.forwardFunction = ghmmwrapper.foba_logp
        self.forwardAlphaFunction = ghmmwrapper.foba_forward      
        self.backwardBetaFunction = ghmmwrapper.foba_backward 
        self.getStatePtr = ghmmwrapper.get_stateptr 
        self.fileWriteFunction = ghmmwrapper.call_model_print
      #  self.getModelPtr = ghmmwrapper.get_model_ptr
      
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
                strout += "\ntransition to state " + str( ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.out_a,i))
            strout +=  "\nIngoing transitions:"
            for i in range(state.in_states):
                strout +=  "\ntransition from state " + str( ghmmwrapper.get_arrayint(state.in_id,i) ) + " with probability " + str(ghmmwrapper.get_arrayd(state.in_a,i))
                strout += "\nint fix:" + str(state.fix) + "\n"
            strout += "Silent states: \n"
            for k in range(hmm.N):
                strout += str(ghmmwrapper.get_arrayint(self.cmodel.silent,k)) + ", "
                strout += "\n"
        return strout
    

    def setEmission(self, i, distributionParemters):
        """ Set the emission distribution parameters for a discrete model."""
        assert len(distributionParemters) == self.N
        state = self.getStatePtr(self.cmodel.s,i)
        for i in range(M):
            ghmmwrapper.set_arrayd(state.b,i,distributionParameters[i])

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
            

class GaussianEmissionHMM(HMM):
    """ GaussianEmissionHMM are HMMs which have a Gaussian distribution
        per state determining the emission probababilities
    """
  
    def __init__(self, emissionDomain, distribution, cmodel):
        HMM.__init__(self, emissionDomain, distribution, cmodel)
        self.silent = 0  # flag indicating whether the model type does include silent states
        
        self.samplingFunction = ghmmwrapper.smodel_generate_sequences
        self.viterbiFunction = ghmmwrapper.sviterbi
        self.forwardFunction = ghmmwrapper.sfoba_logp
        self.forwardAlphaFunction = ghmmwrapper.sfoba_forward
        self.backwardBetaFunction = ghmmwrapper.sfoba_backward 
        self.getStatePtr = ghmmwrapper.get_sstate_ptr
        self.fileWriteFunction = ghmmwrapper.call_smodel_print
        #  self.getModelPtr = ghmmwrapper.get_model_ptr

    def loglikelihood_sqd(self, sequenceSet):
        # XXX REMOVE soon XXX. Currently needed in GLQQuery.py
        """ Compute log( P[emissionSequences| model]) using the forward algorithm

            sequences can either be a SequenceSet or a EmissionSequence

            Result: log( P[emissionSequences| model]) of type float
                    numarray vector of floats
            
            Note: The implementation will not compute the full forward matrix 
        """
        likelihood = ghmmwrapper.double_array(len(sequences))
        result = ghmmwrapper.smodel_individual_likelihoods(self.cmodel, sequences, likelihood)
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
        hmm = self.cmodel
        strout = "\nOverview of HMM:"
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
            strout += "\nOutgoing transitions:"

            for i in range( state.out_states):
                strout += "\ntransition to state " + str(ghmmwrapper.get_arrayint(state.out_id,i) ) + " with probability = " + str(ghmmwrapper.get_2d_arrayd(state.out_a,0,i))
            strout +=  "\nIngoing transitions:"
            for i in range(state.in_states):
                strout += "\ntransition from state " + str(ghmmwrapper.get_arrayint(state.in_id,i) ) +" with probability = "+ str(ghmmwrapper.get_2d_arrayd(state.in_a,0,i))
            strout += "\nint fix:" + str(state.fix) + "\n"
        return strout

    
    # different function signatures require overloading of parent class methods    
    def sample(self, seqNr ,T):
        """ Sample emission sequences 


        """
        seqPtr = self.samplingFunction(self.cmodel,0,T,seqNr,0,-1) # XXX label assignment
        return SequenceSet(self.emissionDomain,seqPtr)
        

    def sampleSingle(self, T):
        """ Sample a single emission sequence of length at most T.
            Returns a Sequence object.
        """
        seqPtr = self.samplingFunction(self.cmodel,0,T,1,0,-1) # XXX label assignment
        return EmissionSequence(self.emissionDomain,seqPtr)
        
    def forward(self, emissionSequence):
        """

            Result: the (N x T)-matrix containing the forward-variables
                    and the scaling vector
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise AttributeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)

        logP = ghmmwrapper.double_array(1)		
        i = self.cmodel.N


        t = ghmmwrapper.get_arrayint(emissionSequence.cseq.seq_len,0)
        calpha = ghmmwrapper.matrix_d_alloc(t,i)
        cscale = ghmmwrapper.double_array(t)

        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)	
		
        error = self.forwardAlphaFunction(self.cmodel, seq,t, None, calpha, cscale, logP)
        if error == -1:
            print "ERROR: Forward finished with -1: Sequence " + str(seq_nr) + " cannot be build."
            exit
        
        # translate alpha / scale to python lists 
        pyscale = ghmmhelper.arrayd2list(cscale, t)
        pyalpha = ghmmhelper.matrixd2list(calpha,t,i)
        
        ghmmwrapper.freearray(logP)
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(calpha,t)
        return (pyalpha,pyscale)         
        
    def backward(self, emissionSequence, scalingVector):
        """

            Result: the (N x T)-matrix containing the backward-variables
        """
        if not isinstance(emissionSequence,EmissionSequence):
            raise AttributeError, "EmissionSequence required, got " + str(emissionSequence.__class__.__name__)
        
        seq = emissionSequence.getPtr(emissionSequence.cseq.seq,0)
        
        # parsing 'scalingVector' to C double array.
        cscale = ghmmhelper.list2arrayd(scalingVector)
        
        # alllocating beta matrix
        t = ghmmwrapper.get_arrayint(emissionSequence.cseq.seq_len,0)
        cbeta = ghmmwrapper.double_2d_array(t, self.cmodel.N)
        
        error = self.backwardBetaFunction(self.cmodel,seq,t,None,cbeta,cscale)
        if error == -1:
            print "ERROR: backward finished with -1: EmissionSequence cannot be build."
            exit
        
        pybeta = ghmmhelper.matrixd2list(cbeta,t,self.cmodel.N)

        # deallocation
        ghmmwrapper.freearray(cscale)
        ghmmwrapper.free_2darrayd(cbeta,t)
        return pybeta        


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
                    if inId == i:
                        ghmmwrapper.set_arrayd(inState.in_a,k,normP) # updating in probabilities



