#!/usr/bin/env python2.3
################################################################################
#
#       This file is part of GHMM (General Hidden Markov Model library) 
#
#       file:    ghmmunittests.py
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
Testing GHMM

- EmissionDomainTests
- AlphabetTests
- FloatTestsTests

- AbstractDistributionTests
- DiscreteDistributionTests
- GaussianDistributionTests
- GaussianMixtureDistributionTests

- EmissionSequenceTests
- SequenceSet

...

"""


import unittest
import ghmm
import random

class AlphabetTests(unittest.TestCase):
    """Unittests for Emissiondomains subclasses"""

    def setUp(self):
        self.binaryAlphabet = ghmm.Alphabet(['zero','one'])
        self.dna = ['a','c','g','t']
        self.dnaAlphabet = ghmm.Alphabet(self.dna)


    def testinternalexternal(self):
        # Check that internal -> external is a bijection
        for l in self.dna:
            self.assertEqual(l, self.dnaAlphabet.external(self.dnaAlphabet.internal(l)))

        self.assertRaises(KeyError, self.dnaAlphabet.internal, '')

        self.assertRaises(KeyError, self.dnaAlphabet.internal, 'x')

        self.assertRaises(KeyError, self.dnaAlphabet.external, -1)
        
        self.assertRaises(KeyError, self.dnaAlphabet.external, len(self.dna) + 1)
    

    def testinternalexternalSequence(self):
        extseq = ['a','c','g','t','a','g','t']
        intseq = self.dnaAlphabet.internalSequence(extseq)
        self.assertEqual(min(intseq), 0)
        self.assertEqual(max(intseq), len(self.dna) - 1)
        result = self.dnaAlphabet.externalSequence(intseq)
        for i in range(len(result)):
            self.assertEqual(extseq[i], result[i])
            self.assertEqual(intseq[i], self.dnaAlphabet.internal(extseq[i]))
            self.assertEqual(result[i], self.dnaAlphabet.external(intseq[i]))

        
        extseq = ['a','c','g','x','a','g','t']
        self.assertRaises(KeyError, self.dnaAlphabet.internalSequence, extseq)
        intseq[3] = 5
        self.assertRaises(KeyError, self.dnaAlphabet.externalSequence, intseq)
       
        
                            
    def testlen(self):
        self.assertEqual(len(self.binaryAlphabet),2)
        self.assertEqual(len(self.dnaAlphabet),len(self.dna))

        self.assertEqual(self.binaryAlphabet.size(),2)
        self.assertEqual(self.dnaAlphabet.size(),len(self.dna))
    

class EmissionSequenceTests(unittest.TestCase):
    
    def setUp(self):
        i_alph = ghmm.IntegerRange(0,5)
        d_alph = ghmm.Float()
        self.i_seq = ghmm.EmissionSequence(i_alph,[1,2,0,0,0,3,4])
        self.d_seq = ghmm.EmissionSequence(d_alph,[1.3, 2.1, 0.8, 0.1, 0.03, 3.6, 43.3])
        
        
    def testprint(self):
        s = "\nEmissionSequence Instance:\nlength 7, weight 1.0:\n1200034"
        self.assertEqual(str(self.i_seq),s)
        
        s2 = "\nEmissionSequence Instance:\nlength 7, weight 1.0:\n1.3 2.1 0.8 0.1 0.03 3.6 43.3 "
        self.assertEqual(str(self.d_seq),s2)
        
        
    def testattributes(self):
        self.assertEqual(self.i_seq.cseq.state_labels,None)
        self.assertEqual(self.i_seq.cseq.state_labels_len,None)    
        self.assertEqual(self.i_seq.cseq.seq_number,1)    
        self.assertEqual(len(self.i_seq),7)    
    
        self.assertEqual(self.d_seq.cseq.seq_number,1)    
        self.assertEqual(len(self.d_seq),7) 
        
    def testitemaccess(self):
        b = self.i_seq[5]
        self.assertEqual(b,3)    
        self.i_seq[5] = 1
        self.assertEqual(self.i_seq[5],1)            
        
        b2 = self.d_seq[1]
        self.assertEqual(b2,2.1)    
        self.d_seq[1] = 8.34
        self.assertEqual(self.d_seq[1],8.34) 
    
    def testwrite(self):
        self.i_seq.write("ghmmunittest_testwrite.seq")
        self.d_seq.write("ghmmunittest_testwrite.seq")
        
    def testweightaccess(self):
        w = self.i_seq.getWeight()
        self.assertEqual(w,1.0)
        self.i_seq.setWeight(4.0)
        w = self.i_seq.getWeight()
        self.assertEqual(w,4.0)
        
        w2 = self.d_seq.getWeight()
        self.assertEqual(w2,1.0)
        self.d_seq.setWeight(2.0)
        w2 = self.d_seq.getWeight()
        self.assertEqual(w2,2.0)
        

class SequenceSetTests(unittest.TestCase):
    
    def setUp(self):
        self.i_alph = ghmm.IntegerRange(0,7)
        self.d_alph = ghmm.Float()
        self.i_seq = ghmm.SequenceSet(self.i_alph,[ [1,2,3,4,5],[0,3,0],[4,3,2,2,1,1,1,1], [0,0,0,2,1],[1,1,1,1,1,1] ])
        self.d_seq = ghmm.SequenceSet(self.d_alph,[ [1.5,2.3,3.7,4.1,5.1],[0.0,3.1,0.7],[4.4,3.05,2.0,2.4,1.2,1.8,1.0,1.0], [0.4,0.1,0.33,2.7,1.345],[1.0,1.0,1.0,1.0,1.0,1.0] ])
 
 
    def testprint(self):
        s = "\nNumber of sequences: 5\nSeq 0, length 5, weight 1.0:\n12345\nSeq 1, length 3, weight 1.0:\n030\nSeq 2, length 8, weight 1.0:\n43221111\nSeq 3, length 5, weight 1.0:\n00021\nSeq 4, length 6, weight 1.0:\n111111"
        self.assertEqual(str(self.i_seq),s)

        s2 = "\nNumber of sequences: 5\nSeq 0, length 5, weight 1.0:\n1.5 2.3 3.7 4.1 5.1 \nSeq 1, length 3, weight 1.0:\n0.0 3.1 0.7 \nSeq 2, length 8, weight 1.0:\n4.4 3.05 2.0 2.4 1.2 1.8 1.0 1.0 \nSeq 3, length 5, weight 1.0:\n0.4 0.1 0.33 2.7 1.345 \nSeq 4, length 6, weight 1.0:\n1.0 1.0 1.0 1.0 1.0 1.0 "
        self.assertEqual(str(self.d_seq),s2)

       
    def testattributes(self):
        self.assertEqual(len(self.i_seq),5)
        self.assertEqual(self.i_seq.sequenceLength(1),3)
        
        self.assertEqual(len(self.d_seq),5)
        self.assertEqual(self.d_seq.sequenceLength(4),6)
     
    def testgetitem(self):
        s = self.i_seq[2]
        self.assertEqual(len(s),8)
        
        s2 = self.d_seq[4]
        self.assertEqual(len(s2),6)
        
    
    def testweightaccess(self):
        w = self.i_seq.getWeight(4)
        self.assertEqual(w,1.0)
        self.i_seq.setWeight(4,4.0)
        w = self.i_seq.getWeight(4)
        self.assertEqual(w,4.0)
        
        w2 = self.d_seq.getWeight(2)
        self.assertEqual(w2,1.0)
        self.d_seq.setWeight(2,7.0)
        w2 = self.d_seq.getWeight(2)
        self.assertEqual(w2,7.0)
        

    def testmerge(self):
        """Merging two SequenceSets   """
        wrong = 4  # wrong argument type to merge
        self.assertRaises(TypeError,self.i_seq.merge,wrong)

        mseq = ghmm.SequenceSet(self.i_alph,[ [1,4,0,4,5,3],[1,2,3,0] ])        
        self.i_seq.merge(mseq)
        self.assertEqual(len(self.i_seq),7)
        s = "\nNumber of sequences: 7\nSeq 0, length 5, weight 1.0:\n12345\nSeq 1, length 3, weight 1.0:\n030\nSeq 2, length 8, weight 1.0:\n43221111\nSeq 3, length 5, weight 1.0:\n00021\nSeq 4, length 6, weight 1.0:\n111111\nSeq 5, length 6, weight 1.0:\n140453\nSeq 6, length 4, weight 1.0:\n1230"
        self.assertEqual(str(self.i_seq),s)

        d_mseq = ghmm.SequenceSet(self.d_alph,[ [7.5,4.0,1.2],[0.4,0.93,3.3,2.54] ])    
        self.d_seq.merge(d_mseq)
        self.assertEqual(len(self.d_seq),7)
        s2 = "\nNumber of sequences: 7\nSeq 0, length 5, weight 1.0:\n1.5 2.3 3.7 4.1 5.1 \nSeq 1, length 3, weight 1.0:\n0.0 3.1 0.7 \nSeq 2, length 8, weight 1.0:\n4.4 3.05 2.0 2.4 1.2 1.8 1.0 1.0 \nSeq 3, length 5, weight 1.0:\n0.4 0.1 0.33 2.7 1.345 \nSeq 4, length 6, weight 1.0:\n1.0 1.0 1.0 1.0 1.0 1.0 \nSeq 5, length 3, weight 1.0:\n7.5 4.0 1.2 \nSeq 6, length 4, weight 1.0:\n0.4 0.93 3.3 2.54 "
        self.assertEqual(str(self.d_seq),s2)        
        
    
    def testgetsubset(self):
        i_subseq = self.i_seq.getSubset([2,1,3])
        s = "\nNumber of sequences: 3\nSeq 0, length 8, weight 1.0:\n43221111\nSeq 1, length 3, weight 1.0:\n030\nSeq 2, length 5, weight 1.0:\n00021"
        self.assertEqual(str(i_subseq),s)
        self.assertEqual(len(i_subseq),3)
        self.assertEqual(i_subseq.sequenceLength(0),8)
        
        d_subseq = self.d_seq.getSubset([0,4])
        s2 = "\nNumber of sequences: 2\nSeq 0, length 5, weight 1.0:\n1.5 2.3 3.7 4.1 5.1 \nSeq 1, length 6, weight 1.0:\n1.0 1.0 1.0 1.0 1.0 1.0 "
        self.assertEqual(str(d_subseq),s2)
        self.assertEqual(len(d_subseq),2)
        self.assertEqual(d_subseq.sequenceLength(0),5)
        
    def testwrite(self):
       self.i_seq.write("ghmmunittest_testwrite.seq") 
       self.d_seq.write("ghmmunittest_testwrite.seq") 

class DiscreteEmissionHMMTests(unittest.TestCase):
    def setUp(self):
        self.model = ghmm.HMMFromMatrices(ghmm.DNA,ghmm.DiscreteDistribution(ghmm.DNA),
                       [[0.3,0.3,0.4],[0.6,0.1,0.3],[1.0,0.0,0.0]],
                       [[0.0,0.5,0.5,0.0],[0.1,0.0,0.8,0.1], [0.0,0.0,0.0,0.0]],
                       [1.0,0,0])
                       
    def testaccessfunctions(self):
        self.assertEqual(self.model.N,3)
        self.assertEqual(self.model.M,4)
        
        pi = self.model.getInitial(2)
        self.assertEqual(pi,0)
        self.model.setInitial(2,0.5,fixProb=1)
        pi = self.model.getInitial(2)
        self.assertEqual(pi,0.5)
        
        trans = self.model.getTransition(0,1)
        self.assertEqual(trans, 0.3)
        self.model.setTransition(0,1,0.6)
        trans = self.model.getTransition(0,1)
        self.assertEqual(trans, 0.6)
        
        emission = self.model.getEmission(1)
        self.assertEqual(emission, [0.1, 0.0, 0.8, 0.1] )
        self.model.setEmission(1,[0.2,0.2,0.2,0.4])
        emission = self.model.getEmission(1)
        self.assertEqual(emission,[0.2,0.2,0.2,0.4] )        

   
class GaussianEmissionHMMTests(unittest.TestCase):

    def setUp(self):
        self.oneState = self.oneStateModel(0.0, 0.2)

    def oneStateModel(self, mean, var):
        # one state model with N(mean, var)
        return ghmm.HMMFromMatrices(ghmm.Float(),
                                    ghmm.GaussianDistribution(ghmm.Float),
                                    [[1.0]],[[mean, var]], [1.0])

    def testdel(self):
        del(self.oneState)


class DiscreteModelHMMTests(unittest.TestCase):

    def setUp(self):
        self.linearmodel = self.oneModel()

    def oneModel(self):
        A  = [[0.3, 0.6,0.1],[0.0, 0.5, 0.5],[0.0,0.0,1.0]]
        B  = [[0.5, 0.5],[0.5,0.5],[1.0,0.0]]
        pi = [1.0, 0.0, 0.0]
        return ghmm.HMMFromMatrices(ghmm.IntegerRange(0,2),
                                    ghmm.DiscreteDistribution(ghmm.IntegerRange(0,2)),
                                    A, B, pi)
        
    def testdel(self):
        """  test for explicit construction and destruction """
        del self.linearmodel
        print "===== testing for construction/destruction ===="
        for i in range(100):
            mo = self.oneModel()

    # TO DO: testing XML-file read
    
if __name__ == '__main__':
    unittest.main()
