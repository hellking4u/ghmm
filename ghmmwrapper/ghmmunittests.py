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
import ghmmwrapper
import random
import ghmmwrapper

class AlphabetTests(unittest.TestCase):
    """Unittests for Emissiondomains subclasses"""

    def setUp(self):
        self.binaryAlphabet = ghmm.Alphabet(['zero','one'])
        self.dna = ['a','c','g','t']
        self.dnaAlphabet = ghmm.Alphabet(self.dna)


    def testinternalexternal(self):
        """ Check that internal -> external is a bijection """
        
        # print"\ntestinternalexternal ",
        for l in self.dna:
            self.assertEqual(l, self.dnaAlphabet.external(self.dnaAlphabet.internal(l)))

        self.assertRaises(KeyError, self.dnaAlphabet.internal, '')

        self.assertRaises(KeyError, self.dnaAlphabet.internal, 'x')

        self.assertRaises(KeyError, self.dnaAlphabet.external, -1)
        
        self.assertRaises(KeyError, self.dnaAlphabet.external, len(self.dna) + 1)
    

    def testinternalexternalSequence(self):
        """ Check internal -> external applied to a sequence """
        # print"\ntestinternalexternalSequence ",
        
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
        # print"\ntestlen ",
        
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
        # print"\ntestprint ",
        s = "\nEmissionSequence Instance:\nlength 7, weight 1.0:\n1200034"
        self.assertEqual(str(self.i_seq),s)
        
        s2 = "\nEmissionSequence Instance:\nlength 7, weight 1.0:\n1.3 2.1 0.8 0.1 0.03 3.6 43.3 "
        self.assertEqual(str(self.d_seq),s2)
        
        
    def testattributes(self):
        # print"\ntestattributes ", 
        self.assertEqual(self.i_seq.cseq.state_labels,None)
        self.assertEqual(self.i_seq.cseq.state_labels_len,None)    
        self.assertEqual(self.i_seq.cseq.seq_number,1)    
        self.assertEqual(len(self.i_seq),7)    
    
        self.assertEqual(self.d_seq.cseq.seq_number,1)    
        self.assertEqual(len(self.d_seq),7) 
        
    def testitemaccess(self):
        # print"\ntestitemaccess ",
        b = self.i_seq[5]
        self.assertEqual(b,3)    
        self.i_seq[5] = 1
        self.assertEqual(self.i_seq[5],1)            
        
        b2 = self.d_seq[1]
        self.assertEqual(b2,2.1)    
        self.d_seq[1] = 8.34
        self.assertEqual(self.d_seq[1],8.34) 
    
    def testwrite(self):
        # print"\ntestwrite ",
        self.i_seq.write("ghmmunittests_testwrite.seq")
        self.d_seq.write("ghmmunittests_testwrite.seq")
        
    def testweightaccess(self):
        # print"\ntestweightaccess ",
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

    def testlabelaccess(self):
        self.i_seq.getLabel()   
        l = self.d_seq.getLabel()
        self.assertEqual(l,-1)
        self.d_seq.setLabel(5)
        l = self.d_seq.getLabel()
        self.assertEqual(l,5)


class SequenceSetTests(unittest.TestCase):
    
    def setUp(self):
        self.i_alph = ghmm.IntegerRange(0,7)
        self.d_alph = ghmm.Float()
        self.i_seq = ghmm.SequenceSet(self.i_alph,[ [1,2,3,4,5],[0,3,0],[4,3,2,2,1,1,1,1], [0,0,0,2,1],[1,1,1,1,1,1] ])
        self.d_seq = ghmm.SequenceSet(self.d_alph,[ [1.5,2.3,3.7,4.1,5.1],[0.0,3.1,0.7],[4.4,3.05,2.0,2.4,1.2,1.8,1.0,1.0], [0.4,0.1,0.33,2.7,1.345],[1.0,1.0,1.0,1.0,1.0,1.0] ])
 
 
    def testprint(self):
        # print"\ntestprint ",
        s = "\nNumber of sequences: 5\nSeq 0, length 5, weight 1.0:\n12345\nSeq 1, length 3, weight 1.0:\n030\nSeq 2, length 8, weight 1.0:\n43221111\nSeq 3, length 5, weight 1.0:\n00021\nSeq 4, length 6, weight 1.0:\n111111"
        self.assertEqual(str(self.i_seq),s)

        s2 = "\nNumber of sequences: 5\nSeq 0, length 5, weight 1.0:\n1.5 2.3 3.7 4.1 5.1 \nSeq 1, length 3, weight 1.0:\n0.0 3.1 0.7 \nSeq 2, length 8, weight 1.0:\n4.4 3.05 2.0 2.4 1.2 1.8 1.0 1.0 \nSeq 3, length 5, weight 1.0:\n0.4 0.1 0.33 2.7 1.345 \nSeq 4, length 6, weight 1.0:\n1.0 1.0 1.0 1.0 1.0 1.0 "
        self.assertEqual(str(self.d_seq),s2)

       
    def testattributes(self):
        # print"\ntestattributes ",
        self.assertEqual(len(self.i_seq),5)
        self.assertEqual(self.i_seq.sequenceLength(1),3)
        
        self.assertEqual(len(self.d_seq),5)
        self.assertEqual(self.d_seq.sequenceLength(4),6)
     
    def testgetitem(self):
        # print"\ntestgetitem ",
        s = self.i_seq[2]
        self.assertEqual(len(s),8)
        
        s2 = self.d_seq[4]
        self.assertEqual(len(s2),6)
        
    
    def testweightaccess(self):
        # print"\ntestweightaccess ",
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
        # print"\ntestmerge ",
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
        # print"\ntestgetsubset ",
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
       # print"\ntestwrite ",
       self.i_seq.write("ghmmunittests_testwrite.seq") 
       self.d_seq.write("ghmmunittests_testwrite.seq") 
       
       
       
    def testlabelaccess(self):
       self.i_seq.getLabel(2)   
       l = self.d_seq.getLabel(3)
       self.assertEqual(l,-1)
       self.d_seq.setLabel(3,8)
       l = self.d_seq.getLabel(3)
       self.assertEqual(l,8)
  

class DiscreteEmissionHMMTests(unittest.TestCase):
    def setUp(self):
        
        self.A = [[0.3,0.3,0.4],[0.6,0.1,0.3],[1.0,0.0,0.0]]
        self.B = [[0.0,0.5,0.5,0.0],[0.1,0.0,0.8,0.1], [0.0,0.0,0.0,0.0]]
        self.pi = [1.0,0,0]
        self.model = ghmm.HMMFromMatrices(ghmm.DNA,ghmm.DiscreteDistribution(ghmm.DNA), self.A, self.B, self.pi)
                       
    def testaccessfunctions(self):

        # print"\ntestaccessfunctions",

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
        
        # introducing silent state
        self.model.setEmission(1,[0.0,0.0,0.0,0.0])
        emission = self.model.getEmission(1)
        self.assertEqual(emission,[0.0,0.0,0.0,0.0] ) 
        self.assertEqual(self.model.cmodel.model_type,4)
        self.assertEqual(ghmmwrapper.get_arrayint(self.model.cmodel.silent,1),1)
        
        # removing silent state
        self.model.setEmission(1,[0.2,0.2,0.2,0.4])
        emission = self.model.getEmission(1)
        self.assertEqual(emission,[0.2,0.2,0.2,0.4] )  
        self.assertEqual(self.model.cmodel.model_type,4)
        self.assertEqual(ghmmwrapper.get_arrayint(self.model.cmodel.silent,1),0)
        
        # removing last silent state
        self.model.setEmission(2,[0.25,0.25,0.25,0.25])
        emission = self.model.getEmission(2)
        self.assertEqual(emission,[0.25,0.25,0.25,0.25])  
        self.assertEqual(self.model.cmodel.model_type,0)
        self.assertEqual(ghmmwrapper.get_arrayint(self.model.cmodel.silent,2),0)
    
        # insertin silent state
        self.model.setEmission(2,[0.0,0.0,0.0,0.0])
        emission = self.model.getEmission(2)
        self.assertEqual(emission,[0.0,0.0,0.0,0.0])  
        self.assertEqual(self.model.cmodel.model_type,4)
        self.assertEqual(ghmmwrapper.get_arrayint(self.model.cmodel.silent,2),1)

    
    def getModel(self):
        A  = [[0.3, 0.6,0.1],[0.0, 0.5, 0.5],[0.0,0.0,1.0]]
        B  = [[0.5, 0.5],[0.5,0.5],[1.0,0.0]]
        pi = [1.0, 0.0, 0.0]
        return ghmm.HMMFromMatrices(ghmm.IntegerRange(0,2),
                                    ghmm.DiscreteDistribution(ghmm.IntegerRange(0,2)),
                                    A, B, pi)
        
    def testdel(self):
        """  test for explicit construction and destruction """
        # print"\ntestdel ",
        del self.model
        # print "===== testing for construction/destruction ===="
        for i in range(100):
            mo = self.getModel()
    
    def testtomatrices(self):
        # print"\ntesttomatrices",
        tA,tB,tpi = self.model.toMatrices()
        
        self.assertEqual(self.A,tA)
        self.assertEqual(self.B,tB)
        self.assertEqual(self.pi,tpi)        
            
    def testsample(self):
        # print"\ntestsample ",
        seq = self.model.sampleSingle(100,seed=3586662)
        seq2 = self.model.sample(10,100,seed=3586662)


    def testbaumwelch(self):
        # print"\ntestbaumwelch ",
        seq = self.model.sample(100,100,seed=3586662)
        self.model.baumWelch(seq,5,0.01)
        
        self.model.setEmission(2,[0.25,0.25,0.25,0.25])
        self.model.cmodel.model_type = 1
        ghmmwrapper.set_arrayint(self.model.cmodel.silent,2,0)
        self.model.baumWelch(seq,5,0.01)
        
    def testviterbi(self):
        # print"\ntestviterbi ",
        seq = self.model.sampleSingle(15,seed=3586662)
        path = self.model.viterbi(seq)
        truePath = [0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 2, 0, 1, 0, 1]
        self.assertEqual(path,truePath)    
    
        seq2 = self.model.sample(10,15,seed=3586662)
        path2 = self.model.viterbi(seq2)
        truePath2 = [[0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 2, 0, 1, 0, 1], [0, 1, 0, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0], [0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0], [0, 2, 0, 2, 0, 2, 0, 1, 0, 1, 0, 2, 0, 1, 1, 0, 2, 0, 1, 0], [0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 1, 0, 1], [0, 2, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 2, 0, 1, 0, 1, 0, 1], [0, 2, 0, 2, 0, 2, 0, 2, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0], [0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2, 0, 2, 0, 1, 0, 1, 0], [0, 2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 2, 0, 1, 0], [0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 2, 0, 1, 0, 1]]
        self.assertEqual(path2,truePath2)    

    def testloglikelihood(self):
        # print"\ntestloglikelihood",
        seq = self.model.sampleSingle(100,seed=3586662)
        logp = self.model.loglikelihood(seq)
        self.assert_(logp-93.1053904716 < 10^-8, "Different results in loglikelihood ")

    def testlogprob(self):
        # print"\ntestlogprob ",
        seq = self.model.sampleSingle(15,seed=3586662)
        path = self.model.viterbi(seq)
        logp = self.model.logprob(seq,path)
        self.assert_(logp - 22.4303246929 < 10^-8, "Different results in logprob ")
        
    def testfoba(self):
        
        # print"\ntestfoba ",
        seq = self.model.sampleSingle(40)
        (alpha,scale) = self.model.forward(seq)
        
        talpha = [[0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.48775986167826385, 0.24395091819263898, 0.26828922012909723], [0.0, 0.76923076923076927, 0.23076923076923078], [0.61307901907356943, 0.10899182561307903, 0.27792915531335155], [0.46113149992850677, 0.27262761546160807, 0.26624088460988515], [0.49843909806377906, 0.23245020208516096, 0.2691106998510599], [0.48925773848834581, 0.24233782008947377, 0.26840444142218045], [0.49151597526863366, 0.23990587278762532, 0.26857815194374107], [0.49096046411164368, 0.24050411557207602, 0.2685354203162803], [0.49109711143616286, 0.24035695691490147, 0.26854593164893559], [0.49106349798068222, 0.2403931560208038, 0.26854334599851404], [0.49107176643407435, 0.24038425153253534, 0.26854398203339036], [0.49106973250582803, 0.2403864419168005, 0.26854382557737139], [0.49107023282473811, 0.24038590311182045, 0.26854386406344138], [0.0, 0.76923076923076927, 0.23076923076923075], [0.61307901907356943, 0.10899182561307903, 0.27792915531335155], [0.46113149992850677, 0.27262761546160807, 0.26624088460988515], [0.49843909806377906, 0.23245020208516096, 0.2691106998510599], [0.48925773848834581, 0.24233782008947377, 0.26840444142218045], [0.49151597526863366, 0.23990587278762532, 0.26857815194374107], [0.0, 0.76923076923076916, 0.23076923076923073], [0.71428571428571419, 0.0, 0.2857142857142857], [0.7142857142857143, 0.0, 0.2857142857142857], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.7142857142857143, 0.0, 0.28571428571428575], [0.7142857142857143, 0.0, 0.28571428571428575], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.7142857142857143, 0.0, 0.2857142857142857], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.7142857142857143, 0.0, 0.2857142857142857], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206]]
        self.assertEqual(alpha,talpha)
        
        tscale = [0.69999999999999996, 0.57285714285714284, 0.56965087281795512, 0.57043689532898478, 0.022193996541956595, 0.56461538461538463, 0.57168937329700276, 0.5699361326914828, 0.5703666049776589, 0.57026066621332705, 0.57028672279156123, 0.57028031304744209, 0.57028188974734029, 0.57028150189977711, 0.5702815973050086, 0.57028157383660572, 0.02227675582061845, 0.56461538461538463, 0.57168937329700276, 0.5699361326914828, 0.5703666049776589, 0.57026066621332705, 0.022287899381715846, 0.48461538461538456, 0.34999999999999998, 0.57285714285714284, 0.40236907730673316, 0.57285714285714284, 0.40236907730673316, 0.34999999999999998, 0.34999999999999998, 0.34999999999999998, 0.57285714285714284, 0.56965087281795512, 0.38953070962658143, 0.57285714285714284, 0.56965087281795512, 0.38953070962658143, 0.57285714285714284, 0.56965087281795512]
        self.assertEqual(scale,tscale)
        
        beta = self.model.backward(seq,scale)
        tbeta = [[4.7023258381970335e-08, 4.9049616680238663e-08, 0.0], [7.6476880486013336e-08, 6.4441989048444662e-08, 6.6750394439163765e-08], [8.8750912271283845e-08, 1.2605202032733066e-07, 7.7899391106213766e-08], [2.2011649936367936e-07, 7.3372166454559791e-08, 1.9293676580713878e-07], [1.6587541968363731e-07, 1.6284216085683637e-07, 0.0], [2.4289833862457837e-07, 2.3842109620511866e-07, 2.1510070859124786e-07], [3.6006136297026362e-07, 3.5355497724012006e-07, 3.149099666605184e-07], [5.3239905136437528e-07, 5.2230051265933343e-07, 4.6706904583339777e-07], [7.8672756150767899e-07, 7.7355627164113795e-07, 6.8966832440908326e-07], [1.1663284814732839e-06, 1.1403771297198445e-06, 1.0226275022771566e-06], [1.7145122143536993e-06, 1.6998533968533431e-06, 1.5032019384574998e-06], [2.5739016260687416e-06, 2.465280493395748e-06, 2.25669865080422e-06], [3.6664866296729516e-06, 3.8244853703637139e-06, 3.2146265694123346e-06], [5.936227041531853e-06, 5.0020643568148531e-06, 5.2046463209454603e-06], [6.896579619294006e-06, 9.7951420679827902e-06, 6.0466440192751379e-06], [1.709996643208889e-05, 5.699988810696297e-06, 1.4992564389769577e-05], [1.2807225714868794e-05, 1.2697725891613876e-05, 0.0], [1.9035783321661597e-05, 1.8232454894197703e-05, 1.6857301306648554e-05], [2.7183129465348012e-05, 2.8354523406566894e-05, 2.3774387573954338e-05], [4.398420608268567e-05, 3.7062569873039537e-05, 3.8586960502901783e-05], [5.110751220623913e-05, 7.2587481104513543e-05, 4.4802335690955294e-05], [0.00012671566938798117, 4.2238556462660399e-05, 0.00011110328740486767], [4.7070434823434923e-05, 9.4140869646869846e-05, 0.0], [0.00015207371250648203, 0.00030414742501296406, 0.00015690144941144974], [0.00035483866251512469, 0.00025908854723326563, 0.00050691237502160675], [0.00032265374981534008, 0.00064530749963068016, 0.00028161798612311482], [0.00086550594401837283, 0.00063195672102928807, 0.001075512499384467], [0.00078700200351784007, 0.0015740040070356801, 0.00068690947937966095], [0.0021111017999601578, 0.0042222035999203156, 0.0026233400117261336], [0.0049259041999070345, 0.0098518083998140691, 0.0070370059998671923], [0.011493776466449746, 0.022987552932899492, 0.016419680666356781], [0.026818811755049406, 0.029554356731528338, 0.038312588221499154], [0.047237233589054815, 0.034490678493595578, 0.04122950562635707], [0.042712271973826509, 0.085424543947653017, 0.037489867927821281], [0.11091827741152124, 0.12223204997326255, 0.14237423991275502], [0.19536557500155205, 0.1426478801598634, 0.17051858167217759], [0.17665106399716793, 0.35330212799433586, 0.15505204365202543], [0.45874009543404987, 0.45169198441535702, 0.58883687999055978], [0.68462986472880094, 0.66707525281267788, 0.59755723854134746], [1.0, 1.0, 0.87773059580615509]]
        self.assertEqual(beta,tbeta)


class EmissionSequenceTests(unittest.TestCase):

    def setUp(self):
        self.model   = self.oneModel()
        self.labeled = ghmm.EmissionSequence(ghmm.DNA, 'acgttgatgga',
                                               ['E','R','T','T','T','E','R','T','T','T','R'],
                                               self.model)
        
    def oneModel(self):
        A  = [[0.3,0.6,0.1], [0.0,0.5,0.5], [0.0,0.0,1.0]]
        B  = [[0.5,0.2,0.2,0.1], [0.7,0.1,0.15,0.05], [1.0,0.0,0.0,0.0]]
        pi = [1.0, 0.0, 0.0]
        labellist = ['E','R','T']
        return ghmm.HMMFromMatrices(ghmm.DNA, ghmm.DiscreteDistribution(ghmm.DNA),
                                    A, B, pi, None, None, labellist)


    def test_errors(self):
        pass


    def test_labeled(self):
        #testing length
        self.assertEqual(len(self.labeled), 11)
        #testing sequence
        sequence = ""
        for i in range(len(self.labeled) ):
            sequence += str( self.labeled.emissionDomain.external(self.labeled[i]) )
        self.assertEqual(sequence,'acgttgatgga')
        label = []
        for i in range(len(self.labeled) ):
            label.append(self.model.int2Label(self.labeled.getSymbol(self.labeled.cseq.state_labels, 0, i)))
        self.assertEqual(label,['E','R','T','T','T','E','R','T','T','T','R'])
        


class SequenceSetTests(unittest.TestCase):

    def setUp(self):
        self.seqList = ['aaaa',
                        'acctttg',
                        'ttgggaaaaaa',
                        'ggggggggggggggtaaatttaa',
                        'gggttccgcggaagggggggggctttta']
        self.labelList = [['E','R','T','T'],
                          ['E','R','T','T','E','R','T'],
                          ['E','R','T','T','R','T','T','R','T','E','T'],
                          ['E','R','T','T','R','T','T','R','T','T','R','T','E','T','R','T','T','R','T','T','R','E','T'],
                          ['E','R','T','T','R','T','T','R','T','T','R','T','E','T','R','T','T','R','T','T','R','T','E','T','R','T','T','R'],]
        self.model   = self.oneModel()
        self.seqSet  = ghmm.SequenceSet(ghmm.DNA, self.seqList, self.labelList, self.model)

    def teardown(self):
        del self.seqList
        del self.labelList
        del self.model
        del self.seqSet
        
    def oneModel(self):
        A  = [[0.3,0.6,0.1], [0.0,0.5, 0.5], [0.0,0.0,1.0]]
        B  = [[0.5,0.2,0.2,0.1], [0.7,0.1,0.15,0.05], [1.0,0.0,0.0,0.0]]
        pi = [1.0, 0.0, 0.0]
        labellist = ['E','R','T']
        return ghmm.HMMFromMatrices(ghmm.DNA, ghmm.DiscreteDistribution(ghmm.DNA),
                                    A, B, pi, None, None, labellist)

    def test_seqSet(self):
        self.assertEqual(len(self.seqSet), 5)

        for i in range(len(self.seqSet)):
            
            # testing length
            self.assertEqual(len(self.seqSet.getSequence(i)), len(self.seqList[i]))

            # testing sequence
            sequence = map(self.seqSet.emissionDomain.external, self.seqSet.getSequence(i))
            seq = ""
            for j in range(len(sequence)):
                seq += sequence[j]
            self.assertEqual(seq, self.seqList[i])

            # testing labels
            label = map(self.model.int2Label, self.seqSet.getLabel(i))
            self.assertEqual(label, self.labelList[i])

            
    def test_seqSError1(self):
        try:
            ghmm.SequenceSet(ghmm.DNA, self.seqList, self.labelList)
        except ghmm.UnknownInputType:
            pass
        else:
            fail("expected a ghmm.UnknownInputType")

    def test_seqSError2(self):
        try:
            ghmm.SequenceSet( ghmm.DNA, 'acccgt', self.seqList, self.oneModel)
        except ghmm.UnknownInputType:
            pass
        else:
            fail("expected a ghmm.UnknownInputType")


class BackwardFunctionTest(unittest.TestCase):

    def setUp(self):
        random.seed(0)
        slength = 45
        self.labels = ['One']*slength
        self.allLabels = ['a','b','c','d','e','f','g']
        #,'h','i','j','k','l','m','n','o']
        self.model  = self.oneModel(['One']*11)
        self.model2 = self.oneModel(self.allLabels)
        
        sequence = ""
        for i in range(slength):
            sequence += random.choice(ghmm.DNA.listOfCharacters)
        self.tSeq  = ghmm.EmissionSequence(ghmm.DNA, sequence, self.labels, self.model)

    def tearDown(self):
        del self.model

    #create a random model with len(LabelList) states and 
    def oneModel(self, LabelList):
        no_states = len(LabelList)
        A = []
        B = []
        pi = []
        pisum = 0
        for i in range(no_states):
            asum = 0
            A_e = []
            #get a random A-row
            for j in range(no_states):
                A_e.append(random.random())
                asum += A_e[-1]
            #normalize this A-row
            for j in range(no_states):
                A_e[j] /= asum
            A.append(A_e)
            
            bsum = 0
            B_e = []
            #get a random B-row
            for j in range(4):
                B_e.append(random.random())
                bsum += B_e[-1]
            #normalize this B-row
            for j in range(4):
                B_e[j] /= bsum
            B.append(B_e)
            
            #get random pi
            pi.append(random.random())
            pisum += pi[-1]

        #normalize pi
        for i in range(no_states):
            pi[i] /= pisum

        return ghmm.HMMFromMatrices(ghmm.DNA, ghmm.DiscreteDistribution(ghmm.DNA),
                                    A, B, pi, None, None, LabelList)
    
    def test_oneLabelCompare(self):
        labelSequence      = self.labels
        (alpha, scale)     = self.model.forward( self.tSeq)
        (b_beta)           = self.model.backward( self.tSeq, scale)
        (bl_logp, bl_beta) = self.model.backwardLabels( self.tSeq, labelSequence, scale)

        #compare beta matrizes from backward and backwardLabels (all states share one label)
        self.assertEqual(b_beta, bl_beta)


    def test_allDifferentLabels(self):
        labelSequence = self.allLabels*4

        sequence = ""
        for i in range(len(labelSequence)):
            sequence += random.choice(ghmm.DNA.listOfCharacters)

        Seq  = ghmm.EmissionSequence(ghmm.DNA, sequence, labelSequence, self.model2)
        
        (fl_logp, alpha, scale) =  self.model2.forwardLabels( Seq, labelSequence)
        (bl_logp, bl_beta)      = self.model2.backwardLabels( Seq, labelSequence, scale)

        #check if the beta matrix is at the appropriated entries 0 or different from 0 
        for i in range(len(bl_beta)):
            i = len(bl_beta)-i-1
            for j in range(len(bl_beta[i])):
                if self.model2.index[labelSequence[i]] == ghmmwrapper.get_stateptr(self.model2.cmodel.s,j).label:
                    self.assertNotEqual(bl_beta[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(bl_beta[i][j]) )
                else:
                    self.assertEqual(bl_beta[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(bl_beta[i][j]))

class ForwardFunctionTest(unittest.TestCase):

    def setUp(self):
        random.seed(0)
        slength = 45
        self.labels = ['One']*slength
        self.allLabels = ['a','b','c','d','e','f','g']
        self.model  = self.oneModel(['One']*11)
        self.model2 = self.oneModel(self.allLabels)
        
        sequence = ""
        for i in range(slength):
            sequence += random.choice(ghmm.DNA.listOfCharacters)
        self.tSeq  = ghmm.EmissionSequence(ghmm.DNA, sequence, self.labels, self.model)

    def tearDown(self):
        del self.model

    #create a random model with len(LabelList) states and 
    def oneModel(self, LabelList):
        no_states = len(LabelList)
        A = []
        B = []
        pi = []
        pisum = 0
        for i in range(no_states):
            asum = 0
            A_e = []
            #get a random A-row
            for j in range(no_states):
                A_e.append(random.random())
                asum += A_e[-1]
            #normalize this A-row
            for j in range(no_states):
                A_e[j] /= asum
            A.append(A_e)
            
            bsum = 0
            B_e = []
            #get a random B-row
            for j in range(4):
                B_e.append(random.random())
                bsum += B_e[-1]
            #normalize this B-row
            for j in range(4):
                B_e[j] /= bsum
            B.append(B_e)
            
            #get random pi
            pi.append(random.random())
            pisum += pi[-1]

        #normalize pi
        for i in range(no_states):
            pi[i] /= pisum

        return ghmm.HMMFromMatrices(ghmm.DNA, ghmm.DiscreteDistribution(ghmm.DNA),
                                    A, B, pi, None, None, LabelList)
    
    def test_oneLabelCompare(self):
        labelSequence          = self.labels
        (alpha, scale)         = self.model.forward(self.tSeq)
        (logp, lalpha, lscale) = self.model.forwardLabels(self.tSeq, labelSequence )

        #compare beta matrizes from backward and backwardLabels (all states share one label)
        self.assertEqual(alpha, lalpha)
        self.assertEqual(scale, lscale)

    def test_allDifferentLabels(self):
        labelSequence = self.allLabels*4

        sequence = ""
        for i in range(len(labelSequence)):
            sequence += random.choice(ghmm.DNA.listOfCharacters)

        Seq  = ghmm.EmissionSequence(ghmm.DNA, sequence, labelSequence, self.model2)
        
        (logp, alpha, scale) =  self.model2.forwardLabels( Seq, labelSequence)

        #check if the beta matrix is at the appropriated entries 0 or different from 0 
        for i in range(len(alpha)):
            i = len(alpha)-i-1
            for j in range(len(alpha[i])):
                if self.model2.index[labelSequence[i]] == ghmmwrapper.get_stateptr(self.model2.cmodel.s,j).label:
                    self.assertNotEqual(alpha[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(alpha[i][j]) )
                else:
                    self.assertEqual(alpha[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(alpha[i][j]))


                
    # TO DO: testing XML-file read
   

class GaussianEmissionHMMTests(unittest.TestCase):

    def setUp(self):
        # print"setUp"
        F = ghmm.Float()
        self.A = [[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]]
        self.B = [[0.0,1.0],[-1.0,0.5], [1.0,0.2]]
        self.pi = [1.0,0.0,0.0]
        self.model = ghmm.HMMFromMatrices(F,ghmm.GaussianDistribution(F), self.A, self.B, self.pi)

    def testaccessfunctions(self):
        # print"\ntestaccessfunctions",
        
        self.assertEqual(self.model.N,3)
        self.assertEqual(self.model.M,1)
        
        pi = self.model.getInitial(2)
        self.assertEqual(pi,0)
        self.model.setInitial(2,0.5,fixProb=1)
        pi = self.model.getInitial(2)
        self.assertEqual(pi,0.5)
        
        trans = self.model.getTransition(0,1)
        self.assertEqual(trans, 1.0)
        self.model.setTransition(0,1,0.6)
        
        trans = self.model.getTransition(0,1)
        self.assertEqual(trans, 0.6)
        
        emission = self.model.getEmission(1)
        self.assertEqual(emission, (-1.0,0.5) )
        self.model.setEmission(1,(3.0,0.5))
        
        emission = self.model.getEmission(1)
        self.assertEqual(emission, (3.0,0.5))
    
    def testtomatrices(self):
        # print"\ntesttomatrices ",
        tA,tB,tpi = self.model.toMatrices()
        
        self.assertEqual(self.A,tA)
        self.assertEqual(self.B,tB)
        self.assertEqual(self.pi,tpi)    
    
    def testsample(self):
        # print"\ntestsample ",
        seq = self.model.sampleSingle(100,seed=3586662)
        seq2 = self.model.sample(10,100,seed=3586662)


    def testbaumwelch(self):
        # print"\ntestbaumwelch",
        seq = self.model.sample(100,100,seed=0)
        self.model.baumWelch(seq,5,0.01)
    
    
    def oneStateModel(self, mean, var):
        # one state model with N(mean, var)
        return ghmm.HMMFromMatrices(ghmm.Float(),
                                    ghmm.GaussianDistribution(ghmm.Float),
                                    [[1.0]],[[mean, var]], [1.0])

    def testdel(self):
        # print"\ntestdel ",
        del(self.model)


if __name__ == '__main__':
    unittest.main()


suite = unittest.TestSuite()
#suite.addTest(AlphabetTests("testinternalexternal"))
#suite.addTest(AlphabetTests("testinternalexternalSequence"))
#suite.addTest(AlphabetTests("testlen"))

#runner = unittest.TextTestRunner()
# runner.run(suite)




