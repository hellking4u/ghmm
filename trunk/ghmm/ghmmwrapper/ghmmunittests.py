#!/usr/bin/env python2.3
################################################################################
#
#       This file is part of the General Hidden Markov Model Library,
#       GHMM version __VERSION__, see http://ghmm.org
#
#       file:    ghmmunittests.py
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
from Numeric import array, sum

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

        self.assertEqual(len(self.binaryAlphabet),2)
        self.assertEqual(len(self.dnaAlphabet),len(self.dna))
    

class EmissionSequenceTests(unittest.TestCase):
    
    def setUp(self):
        i_alph = ghmm.IntegerRange(0,5)
        d_alph = ghmm.Float()
        l_domain = ghmm.LabelDomain(['E','R','T'])
        self.i_seq = ghmm.EmissionSequence(i_alph,[1,2,0,0,0,3,4])
        self.d_seq = ghmm.EmissionSequence(d_alph,[1.3, 2.1, 0.8, 0.1, 0.03, 3.6, 43.3])
        self.labeled = ghmm.EmissionSequence(ghmm.DNA, list('acgttgatgga'),labelDomain=l_domain,
                                            labelInput= ['E','R','T','T','T','E','R','T','T','T','R'])
        
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
        self.i_seq.getSeqLabel()   
        l = self.d_seq.getSeqLabel()
        self.assertEqual(l,-1)
        self.d_seq.setSeqLabel(5)
        l = self.d_seq.getSeqLabel()
        self.assertEqual(l,5)

    def testerrors(self):
        pass


    def testlabeled(self):
        #testing length
        self.assertEqual(len(self.labeled), 11)
        #testing sequence
        sequence = ""
        for i in range(len(self.labeled) ):
            sequence += str( self.labeled.emissionDomain.external(self.labeled[i]) )
        self.assertEqual(sequence,'acgttgatgga')
        label = []
        for i in range(len(self.labeled) ):
            label.append(self.labeled.labelDomain.external(self.labeled.getSymbol(self.labeled.cseq.state_labels, 0, i)))
        self.assertEqual(label,['E','R','T','T','T','E','R','T','T','T','R'])


class SequenceSetTests(unittest.TestCase):
    
    def setUp(self):
        self.i_alph = ghmm.IntegerRange(0,7)
        self.d_alph = ghmm.Float()
        self.l_domain = ghmm.LabelDomain(['E','R','T'])
        
        self.i_seq = ghmm.SequenceSet(self.i_alph,[ [1,2,3,4,5],[0,3,0],[4,3,2,2,1,1,1,1], [0,0,0,2,1],[1,1,1,1,1,1] ])
        self.d_seq = ghmm.SequenceSet(self.d_alph,[ [1.5,2.3,3.7,4.1,5.1],[0.0,3.1,0.7],[4.4,3.05,2.0,2.4,1.2,1.8,1.0,1.0], [0.4,0.1,0.33,2.7,1.345],[1.0,1.0,1.0,1.0,1.0,1.0] ])

        self.seqList = [list('aaaa'),
                        list('acctttg'),
                        list('ttgggaaaaaa'),
                        list('ggggggggggggggtaaatttaa'),
                        list('gggttccgcggaagggggggggctttta')]
        self.labelList = [['E','R','T','T'],
                          ['E','R','T','T','E','R','T'],
                          ['E','R','T','T','R','T','T','R','T','E','T'],
                          ['E','R','T','T','R','T','T','R','T','T','R','T','E','T','R','T','T','R','T','T','R','E','T'],
                          ['E','R','T','T','R','T','T','R','T','T','R','T','E','T','R','T','T','R','T','T','R','T','E','T','R','T','T','R'],]        
 
        self.l_seq  = ghmm.SequenceSet(ghmm.DNA, self.seqList,labelDomain=self.l_domain,labelInput= self.labelList)


    def testlabelseqset(self):
        self.assertEqual(len(self.l_seq), 5)

        for i in range(len(self.l_seq)):
            
            # testing length
            self.assertEqual(len(self.l_seq.getSequence(i)), len(self.seqList[i]))

            # testing sequence
            sequence = map(self.l_seq.emissionDomain.external, self.l_seq.getSequence(i))
            seq = []
            for j in range(len(sequence)):
                seq.append(sequence[j])
            self.assertEqual(seq, self.seqList[i])

            # testing labels
            label =  self.l_seq.getStateLabel(i)
            self.assertEqual(label, self.labelList[i])

    # XXX check different input types
    def testseqerror(self):
        
        # self.assertRaises(ghmm.UnknownInputType,ghmm.SequenceSet,)        
        pass
        
        

    
    def testprint(self):
        # print"\ntestprint ",
        s = "\nNumber of sequences: 5\nSeq 0, length 5, weight 1.0:\n12345\nSeq 1, length 3, weight 1.0:\n030\nSeq 2, length 8, weight 1.0:\n43221111\nSeq 3, length 5, weight 1.0:\n00021\nSeq 4, length 6, weight 1.0:\n111111"
        self.assertEqual(str(self.i_seq),s)

        s2 = "\nNumber of sequences: 5\nSeq 0, length 5, weight 1.0:\n1.5 2.3 3.7 4.1 5.1 \nSeq 1, length 3, weight 1.0:\n0.0 3.1 0.7 \nSeq 2, length 8, weight 1.0:\n4.4 3.05 2.0 2.4 1.2 1.8 1.0 1.0 \nSeq 3, length 5, weight 1.0:\n0.4 0.1 0.33 2.7 1.345 \nSeq 4, length 6, weight 1.0:\n1.0 1.0 1.0 1.0 1.0 1.0 "
        self.assertEqual(str(self.d_seq),s2)

        # XXX str(self.l_seq)

       
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
       self.i_seq.getSeqLabel(2)   
       l = self.d_seq.getSeqLabel(3)
       self.assertEqual(l,-1)
       self.d_seq.setSeqLabel(3,8)
       l = self.d_seq.getSeqLabel(3)
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
        self.assertEqual(self.model.getSilentFlag(1),1)
        
        
        # removing silent state
        self.model.setEmission(1,[0.2,0.2,0.2,0.4])
        emission = self.model.getEmission(1)
        self.assertEqual(emission,[0.2,0.2,0.2,0.4] )  
        self.assertEqual(self.model.cmodel.model_type,4)
        self.assertEqual(self.model.getSilentFlag(1),0)
        
        # removing last silent state
        self.model.setEmission(2,[0.25,0.25,0.25,0.25])
        emission = self.model.getEmission(2)
        self.assertEqual(emission,[0.25,0.25,0.25,0.25])  
        self.assertEqual(self.model.cmodel.model_type,0)
        self.assertEqual(self.model.getSilentFlag(2),0)
    
        # inserting silent state
        self.model.setEmission(2,[0.0,0.0,0.0,0.0])
        emission = self.model.getEmission(2)
        self.assertEqual(emission,[0.0,0.0,0.0,0.0])  
        self.assertEqual(self.model.cmodel.model_type,4)
        self.assertEqual(self.model.getSilentFlag(2),1)

    
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
        #print"\ntestsample ",
        seq = self.model.sampleSingle(100,seed=3586662)
        #print "*************" 
        #print seq
        
        seq2 = self.model.sample(10,100,seed=3586662)


    def testbaumwelch(self):
        print"\n**** testbaumwelch ",
        seq = self.model.sample(100,100,seed=3586662)
        self.model.baumWelch(seq,5,0.01)
        
        self.model.setEmission(2,[0.25,0.25,0.25,0.25])
        self.model.cmodel.model_type = 0
        ghmmwrapper.set_arrayint(self.model.cmodel.silent,2,0)
        self.model.baumWelch(seq,5,0.01)
        
    def testviterbi(self):
        # print"\ntestviterbi ",
        seq = self.model.sampleSingle(15,seed=3586662)
        result = self.model.viterbi(seq)
        trueResult = ([0, 2, 0, 1, 0, 2, 0, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0], -20.286344630126756)
        self.assertEqual(result,trueResult)    
    
        seq2 = self.model.sample(10,15,seed=3586662)
        path2 = self.model.viterbi(seq2)
        truePath2 = ([[0, 2, 0, 1, 0, 2, 0, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0], [0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2, 0], [0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0, 1, 0], [0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 1], [0, 1, 0, 2, 0, 2, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1], [0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 1, 0, 2, 0], [0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 1, 0, 1, 0], [0, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1], [0, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1], [0, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 2, 0, 1, 0, 1, 0]], [-20.286344630126756, -22.953572836708705, -20.286344630126752, -20.691809738234916, -22.183464615012635, -22.953572836708702, -21.777999506904472, -19.516236408430682, -20.104023073332801, -22.365786171806587])
        self.assertEqual(path2,truePath2)    

    def testloglikelihood(self):
        # print"\ntestloglikelihood",
        seq = self.model.sampleSingle(100,seed=3586662)
        logp = self.model.loglikelihood(seq)
        self.assert_(logp-93.1053904716 < 10^-8, "Different results in loglikelihood ")

    def testlogprob(self):
        # print"\ntestlogprob ",
        seq = self.model.sampleSingle(15,seed=3586662)
        path,logp = self.model.viterbi(seq)
        logp = self.model.logprob(seq,path)
        self.assert_(logp - 22.4303246929 < 10^-8, "Different results in logprob ")
        
    def testfoba(self):
        
        # print"\ntestfoba ",
        seq = self.model.sampleSingle(40)
        (alpha,scale) = self.model.forward(seq)
        
        talpha = [[0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.7142857142857143, 0.0, 0.2857142857142857], [0.0, 0.76923076923076916, 0.23076923076923075], [0.61307901907356943, 0.10899182561307905, 0.27792915531335155], [0.46113149992850677, 0.27262761546160807, 0.26624088460988515], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.48775986167826385, 0.24395091819263898, 0.26828922012909723], [0.71428571428571419, 0.0, 0.2857142857142857], [0.43640897755610969, 0.29925187032418948, 0.26433915211970072], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.48775986167826385, 0.24395091819263898, 0.26828922012909723], [0.71428571428571419, 0.0, 0.2857142857142857], [0.7142857142857143, 0.0, 0.2857142857142857], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.7142857142857143, 0.0, 0.2857142857142857], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.7142857142857143, 0.0, 0.28571428571428575], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.7142857142857143, 0.0, 0.2857142857142857], [0.0, 0.76923076923076916, 0.23076923076923075], [0.71428571428571419, 0.0, 0.2857142857142857], [0.7142857142857143, 0.0, 0.2857142857142857]]
        self.assertEqual(alpha,talpha)
        
        tscale =   [0.69999999999999996, 0.57285714285714284, 0.56965087281795512, 0.38953070962658143, 0.027857142857142858, 0.56461538461538452, 0.57168937329700276, 0.39770983270578136, 0.57285714285714284, 0.56965087281795512, 0.57043689532898478, 0.39269141068371188, 0.57285714285714284, 0.56965087281795501, 0.57043689532898478, 0.39269141068371188, 0.34999999999999998, 0.34999999999999998, 0.57285714285714284, 0.40236907730673316, 0.57285714285714284, 0.56965087281795512, 0.38953070962658143, 0.34999999999999998, 0.57285714285714284, 0.40236907730673316, 0.57285714285714284, 0.40236907730673316, 0.34999999999999998, 0.34999999999999998, 0.57285714285714284, 0.40236907730673316, 0.57285714285714284, 0.40236907730673316, 0.57285714285714284, 0.56965087281795512, 0.38953070962658143, 0.027857142857142858, 0.48461538461538456, 0.34999999999999998]
        self.assertEqual(scale,tscale)
        
        beta = self.model.backward(seq,scale)
        tbeta = [[5.8361743879278295e-09, 6.4314698721039723e-09, 0.0], [1.0279528241108607e-08, 7.5056872871586669e-09, 8.9721568189227239e-09], [9.2948289439751381e-09, 1.8589657887950276e-08, 8.1583557469115932e-09], [2.4137475429362165e-08, 8.0458251431207221e-09, 3.0982763146583797e-08], [2.0338793275082881e-08, 2.2413370041550582e-08, 0.0], [3.5308273402649789e-08, 2.5780644071776037e-08, 3.1267544566379247e-08], [3.2040261418666808e-08, 6.4080522837333616e-08, 2.8022439208452213e-08], [8.4951513391116509e-08, 7.9556956751669797e-08, 1.0680087139555602e-07], [1.1741230876651627e-07, 1.2938847886519948e-07, 1.0247957124259525e-07], [2.0564634078137072e-07, 1.5015447104671513e-07, 1.8050208521938813e-07], [1.8620358757316112e-07, 3.7240717514632225e-07, 1.6321138157251645e-07], [4.8747032985648482e-07, 4.5651518615759018e-07, 6.2067862524387043e-07], [6.7373746033358206e-07, 7.4245933891323977e-07, 5.8805015241085719e-07], [1.1800435986695057e-06, 8.6161913553646461e-07, 1.0357603709374248e-06], [1.0684768361553578e-06, 2.1369536723107155e-06, 9.3654253862659187e-07], [2.7972111738181117e-06, 5.5944223476362234e-06, 3.5615894538511928e-06], [6.5268260722422604e-06, 1.3053652144484521e-05, 9.3240372460603729e-06], [1.5229260835231941e-05, 1.1119777752709037e-05, 2.1756086907474201e-05], [1.3847921983963737e-05, 2.7695843967927474e-05, 1.2086714948596779e-05], [3.7146503942020763e-05, 4.0935483602285449e-05, 4.6159739946545791e-05], [6.5427883224379181e-05, 4.7772740132086386e-05, 5.7106631243223726e-05], [5.9160398072061651e-05, 0.0001183207961441233, 5.1926891447919982e-05], [0.00015363194561867477, 0.00030726389123734955, 0.00019720132690687219], [0.00035847453977690781, 0.00026174331475774221, 0.00051210648539558258], [0.00032595984229147394, 0.00065191968458294788, 0.00028450360299754585], [0.00087437440654579088, 0.00063843210636676791, 0.0010865328076382465], [0.00079506607035116129, 0.0015901321407023226, 0.00069394794170300864], [0.0021327333408339131, 0.0042654666816678262, 0.0026502202345038712], [0.0049763777952791305, 0.0099527555905582609, 0.007109111136113044], [0.011611548188984636, 0.0084782732808459244, 0.016587925984263768], [0.010558346539190111, 0.021116693078380221, 0.0092155144357020918], [0.028322347699057757, 0.020679809431058044, 0.035194488463967034], [0.025753427272839365, 0.051506854545678729, 0.022478053729410918], [0.069082551795056221, 0.076129039508656404, 0.08584475757613122], [0.12167834525542114, 0.088844506059513842, 0.10620304448727531], [0.11002250091713514, 0.22004500183427028, 0.09657011528208026], [0.2857142857142857, 0.095238095238095247, 0.3667416697237838], [0.1326530612244898, 0.26530612244897961, 0.0], [0.42857142857142855, 0.8571428571428571, 0.44217687074829937], [1.0, 1.0, 1.4285714285714286]]
        self.assertEqual(beta,tbeta)

    def testtiedstated(self):
        f = lambda x: round(x,15)
        e1 = map(f,self.model.getEmission(0))
        t = (-1,1,1)
        self.model.setTieGroups(t)

        self.model.updateTieGroups()
        em2 = map(f,self.model.getEmission(2))
        self.assertEqual(em2, [0.1, 0.0, 0.8, 0.1])
        
        self.model.setEmission(2,[0.2,0.2,0.2,0.4])
        self.model.updateTieGroups()
        em0 = map(f,self.model.getEmission(0))
        self.assertEqual(em0, [0.0,0.5,0.5,0.0])
        em2 = map(f,self.model.getEmission(2))
        self.assertEqual(em2, [0.15, 0.1, 0.5, 0.25])
        
class BackgroundDistributionTests(unittest.TestCase):
    " Tests for background distributions "

    def setUp(self):
        self.sigma = ghmm.Alphabet(['rot','blau','gruen','gelb'])

        self.model = ghmm.HMMFromMatrices(self.sigma,ghmm.DiscreteDistribution(self.sigma),
                       [[0.3,0.3,0.4],[0.6,0.1,0.3],[1.0,0.0,0.0]],
                       [[0.0,0.5,0.5,0.0],[0.1,0.0,0.8,0.1],
                       [0.25,0.25,0.25,0.25, 0.0,0.5,0.5,0.0, 0.1,0.0,0.8,0.1, 0.1,0.35,0.3,0.25 ]],
                       [1.0,0,0])

        self.bg = ghmm.BackgroundDistribution(self.sigma,[[0.2,0.3,0.1,0.4],
                       [0.1,0.2,0.4,0.3, 0.2,0.3,0.1,0.4, 0.25,0.25,0.25,0.25,0.0,0.5,0.5,0.0 ]]
                  )

    def testprint(self):
        s = str(self.bg)
        self.assertEqual(s,"BackgroundDistribution instance:\nNumber of distributions: 2\n\nGHMM Alphabet:\nNumber of symbols: 4\nExternal: ['rot', 'blau', 'gruen', 'gelb']\nInternal: [0, 1, 2, 3]\n\nDistributions:\n  Order: 0\n  1: [0.20000000000000001, 0.29999999999999999, 0.10000000000000001, 0.40000000000000002]\n  Order: 1\n  2: [0.10000000000000001, 0.20000000000000001, 0.40000000000000002, 0.29999999999999999]\n")

    def testmodelbackgroundaccessfunctions(self):
        self.model.setBackground(self.bg,[0,-1,1])
        # deleting background
        del(self.bg)
        s = str(self.model.background)
        self.assertEqual(s,"BackgroundDistribution instance:\nNumber of distributions: 2\n\nGHMM Alphabet:\nNumber of symbols: 4\nExternal: ['rot', 'blau', 'gruen', 'gelb']\nInternal: [0, 1, 2, 3]\n\nDistributions:\n  Order: 0\n  1: [0.20000000000000001, 0.29999999999999999, 0.10000000000000001, 0.40000000000000002]\n  Order: 1\n  2: [0.10000000000000001, 0.20000000000000001, 0.40000000000000002, 0.29999999999999999]\n")

    def testapplybackground(self):
        self.model.setBackground(self.bg,[0,-1,1])
        self.model.applyBackground([0.1,0.2,0.3])
        print self.model
        
        f = lambda x: round(x,15)
        e1 = map(f,self.model.getEmission(0))
        e2 = map(f,self.model.getEmission(1))
        e3 = map(f,self.model.getEmission(2))
                
        self.assertEqual(e1, [0.02, 0.48, 0.46, 0.04])
        self.assertEqual(e2,[ 0.1, 0.0, 0.8, 0.1])
        self.assertEqual(e3, [ 0.205, 0.235, 0.295, 0.265, 0.06, 0.44, 0.38, 0.12, 0.145, 0.075, 0.635, 0.145, 0.07, 0.395, 0.36, 0.175])
        
        
        

    def testbackgroundtraining(self):
        # XXX test for background distributions
         self.model.setEmission(2,[0.25,0.25,0.25,0.25])
        
    # XXX ...        



class StateLabelHMMTests(unittest.TestCase):

    def setUp(self):
        random.seed(0)
        slength = 45
        self.labels = ['One']*slength
        self.allLabels = ['a','b','c','d','e','f','g']
        self.l_domain= ghmm.LabelDomain(['One','a','b','c','d','e','f','g'])
        

        self.A = [[0.0,0.5,0.5],[0.4,0.2,0.4],[0.3,0.3,0.4]]
        self.B = [[0.2,0.1,0.1,0.6],[0.3,0.1,0.1,0.5],
                  [0.25,0.25,0.25,0.25,   0.0, 0.0, 1.0, 0.0,   0.25,0.25,0.25,0.25,  0.25,0.25,0.25,0.25]]
        self.pi = [1.0,0,0.0]

        self.l_domain2 = ghmm.LabelDomain(['fst','scd','thr'])
        self.model = ghmm.HMMFromMatrices(ghmm.DNA,ghmm.DiscreteDistribution(ghmm.DNA), self.A, self.B, self.pi,labelDomain=self.l_domain2,labelList=['fst','scd','thr'])

        sequence = []
        for i in range(slength):
            sequence.append(random.choice(ghmm.DNA.listOfCharacters))
        self.tSeq  = ghmm.EmissionSequence(ghmm.DNA, sequence, labelDomain=self.l_domain,labelInput=self.labels)


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
                                    A, B, pi, None, self.l_domain, LabelList)
            
    def testsample(self):
        # print"\ntestsample ",
        seq = self.model.sampleSingle(100,seed=3586662)
        seq2 = self.model.sample(10,100,seed=3586662)
        
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
        self.assertEqual(trans, 0.5)
        self.model.setTransition(0,1,0.6)
        trans = self.model.getTransition(0,1)
        self.assertEqual(trans, 0.6)
        
        emission = self.model.getEmission(1)
        self.assertEqual(emission, [0.3,0.1,0.1,0.5] )
        
        # introducing silent state
        self.model.setEmission(1,[0.0,0.0,0.0,0.0])
        emission = self.model.getEmission(1)
        self.assertEqual(emission,[0.0,0.0,0.0,0.0] ) 
        self.assertEqual(self.model.cmodel.model_type & 4, 4)
        self.assertEqual(ghmmwrapper.get_arrayint(self.model.cmodel.silent,1),1)
        
        # removing silent state
        self.model.setEmission(1,[0.2,0.2,0.2,0.4])
        emission = self.model.getEmission(1)
        self.assertEqual(emission,[0.2,0.2,0.2,0.4] )  
        print "model_type = ",self.model.cmodel.model_type
        self.assertEqual(self.model.cmodel.model_type & 4,0)
        self.assertEqual(ghmmwrapper.get_arrayint(self.model.cmodel.silent,1),0)
        
        # inserting silent state
        self.model.setEmission(0,[0.0,0.0,0.0,0.0])
        emission = self.model.getEmission(0)
        self.assertEqual(emission,[0.0,0.0,0.0,0.0])  
        self.assertEqual(self.model.cmodel.model_type & 4,4)
        self.assertEqual(ghmmwrapper.get_arrayint(self.model.cmodel.silent,0),1)

        # label access
        labels = self.model.getLabels()
        self.assertEqual(labels,['fst','scd','thr'])
        self.model.setLabels(['fst','thr','fst'])
        labels = self.model.getLabels()
        self.assertEqual(labels,['fst','thr','fst']) 
    
    def testonelabelcomparebackward(self):
        model  = self.oneModel(['One']*11)
        
        labelSequence      = self.labels
        (alpha, scale)     = model.forward( self.tSeq)
        (b_beta)           = model.backward( self.tSeq, scale)
        (bl_logp, bl_beta) = model.backwardLabels( self.tSeq, labelSequence, scale)

        #compare beta matrizes from backward and backwardLabels (all states share one label)
        self.assertEqual(b_beta, bl_beta)


    def testalldifferentlabelsbackward(self):
        model2 = self.oneModel(self.allLabels)
        
        labelSequence = self.allLabels*4

        sequence = []
        for i in range(len(labelSequence)):
            sequence.append(random.choice(ghmm.DNA.listOfCharacters))

        Seq  = ghmm.EmissionSequence(ghmm.DNA, sequence, self.l_domain,labelSequence)
        
        (fl_logp, alpha, scale) =  model2.forwardLabels( Seq, labelSequence)
        (bl_logp, bl_beta)      = model2.backwardLabels( Seq, labelSequence, scale)

        #check if the beta matrix is at the appropriated entries 0 or different from 0 
        for i in range(len(bl_beta)):
            i = len(bl_beta)-i-1
            for j in range(len(bl_beta[i])):
                if model2.labelDomain.internal(labelSequence[i]) == ghmmwrapper.get_stateptr(model2.cmodel.s,j).label:
                    self.assertNotEqual(bl_beta[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(bl_beta[i][j]) )
                else:
                    self.assertEqual(bl_beta[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(bl_beta[i][j]))

    def testonelabelcompareforward(self):
        model  = self.oneModel(['One']*11)
        
        labelSequence          = self.labels
        (alpha, scale)         = model.forward(self.tSeq)
        (logp, lalpha, lscale) = model.forwardLabels(self.tSeq, labelSequence )

        # compare beta matrizes from backward and backwardLabels (all states share one label)
        # XXX due to rounding errors in the Python floating point representation
        # we have to round for 15 decimal positions
        
        f = lambda x: round(x,12) #XXX
        for i in range(len(alpha)):
            alpha[i] = map(f, alpha[i])
            lalpha[i] = map(f, lalpha[i])        
        
        self.assertEqual(alpha, lalpha)

        scale = map(f, scale)
        lscale = map(f, lscale)        
        self.assertEqual(scale, lscale)

    def testalldifferentlabelsforward(self):
        model2  = self.oneModel(['One']*11)
        
        labelSequence = self.allLabels*4

        sequence = []
        for i in range(len(labelSequence)):
            sequence.append(random.choice(ghmm.DNA.listOfCharacters))

        Seq  = ghmm.EmissionSequence(ghmm.DNA, sequence, self.l_domain,labelSequence)
        
        (logp, alpha, scale) =  model2.forwardLabels( Seq, labelSequence)

        #check if the beta matrix is 0 or different from 0 at the appropriate entries 
        for i in range(len(alpha)):
            i = len(alpha)-i-1
            for j in range(len(alpha[i])):
                if model2.labelDomain.internal(labelSequence[i]) == ghmmwrapper.get_stateptr(model2.cmodel.s,j).label:
                    self.assertNotEqual(alpha[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(alpha[i][j]) )
                else:
                    self.assertEqual(alpha[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(alpha[i][j]))



    def testkbest(self): 
        seq = self.model.sampleSingle(20, seed=3586662)
        print seq
        path = self.model.kbest(seq)
        self.assertEqual(path,(['fst', 'thr', 'thr', 'scd', 'fst', 'scd', 'fst', 'scd', 'thr', 'thr', 'fst', 'thr', 'thr', 'thr', 'thr', 'thr', 'scd', 'fst', 'scd', 'fst'], -35.735009627142446)) 


    def testgradientdescent(self):
        A2 = [[0.3,0.2,0.5],[0.1,0.8,0.1],[0.1,0.4,0.5]]
        B2 = [[0.4,0.2,0.2,0.2],[0.4,0.2,0.2,0.2],
             [0.2,0.1,0.1,0.6,   0.25,0.25,0.25,0.25,   0.5,0.1,0.3,0.1, 0.2,0.1,0.1,0.6]]
        pi2 = [0.5,0.5,0.0] 
    
        model2 = ghmm.HMMFromMatrices(ghmm.DNA,ghmm.DiscreteDistribution(ghmm.DNA), A2, B2, pi2,labelDomain=self.l_domain2,labelList=['fst','scd','thr'])
        
        train = self.model.sample(10,300,seed=3586662)
        model2.gradientSearch(train)

    def testbaumwelch(self):
        # print"\ntestbaumwelch ",
        seq = self.model.sample(100,100,seed=3586662)
        self.model.baumWelchLabels(seq,5,0.01)
        
        self.model.setEmission(2,[0.25,0.25,0.25,0.25])
        self.model.cmodel.model_type = 1
        ghmmwrapper.set_arrayint(self.model.cmodel.silent,2,0)
        self.model.baumWelchLabels(seq,5,0.01)   
    
   
    def testviterbilabels(self):
        seq = self.model.sampleSingle(20,seed=3586662)
        p = self.model.viterbiLabels(seq)
        self.assertEqual(p, (['fst', 'thr', 'thr', 'scd', 'fst', 'scd', 'fst', 'scd', 'thr', 'thr', 'fst', 'thr', 'thr', 'thr', 'thr', 'thr', 'scd', 'fst', 'scd', 'fst'], -39.893892710502115))


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


class GaussianMixtureHMMTests(unittest.TestCase):
    def setUp(self):
        # print"setUp"
        F = ghmm.Float()
        self.A = [[0.25,0.5,0.25],[0.3,0.2,0.5],[0.3,0.3,0.4]]

        self.B = [[ [0.0,1.0,2.0],[1.0,2.5,5.5], [0.5,0.3,0.2]],
             [ [2.0,6.0,1.0],[1.0,0.5,0.7], [0.1,0.5,0.4]],
             [ [4.0,5.0,1.0],[1.0,2.5,2.0], [0.3,0.3,0.4]] ]

        self.pi = [1.0,0.0,0.0]
        
        self.model = ghmm.HMMFromMatrices(F,ghmm.GaussianMixtureDistribution(F), self.A, self.B, self.pi)

    def testcomponentfixing(self):
        f = self.model.getMixtureFix(0)
        self.assertEqual(f,[0,0,0])
        self.model.setMixtureFix(0,[0,1,0])    
        f = self.model.getMixtureFix(0)
        self.assertEqual(f,[0,1,0])
        self.model.setMixtureFix(1,[1,1,1])    
        f = self.model.getMixtureFix(1) 
        self.assertEqual(f,[1,1,1])
        
        # XXX check mu,v,u

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
        

class XMLIOTests(unittest.TestCase):

    def setUp(self):        
        self.A = [[0.3,0.3,0.4],[0.6,0.1,0.3],[1.0,0.0,0.0]]
        self.B = [[0.0,0.5,0.5,0.0],[0.1,0.0,0.8,0.1], [0.0,0.0,0.0,0.0]]
        self.pi = [1.0,0,0]
        self.model = ghmm.HMMFromMatrices(ghmm.DNA,ghmm.DiscreteDistribution(ghmm.DNA), self.A, self.B, self.pi)

        # model with labels
        random.seed(0)
        slength = 45
        self.labels = ['One']*slength
        self.allLabels = ['a','b','c','d','e','f','g']
        self.l_domain= ghmm.LabelDomain(['One','a','b','c','d','e','f','g'])
        

        self.A = [[0.0,0.5,0.5],[0.4,0.2,0.4],[0.3,0.3,0.4]]
        self.B = [[0.2,0.1,0.1,0.6],[0.3,0.1,0.1,0.5],
                  [0.25,0.25,0.25,0.25,   0.0, 0.0, 1.0, 0.0,   0.25,0.25,0.25,0.25,  0.25,0.25,0.25,0.25]]
        self.pi = [1.0,0,0.0]

        self.l_domain2 = ghmm.LabelDomain(['fst','scd','thr'])
        self.label_model = ghmm.HMMFromMatrices(ghmm.DNA,ghmm.DiscreteDistribution(ghmm.DNA), self.A, self.B, self.pi,labelDomain=self.l_domain2,labelList=['fst','scd','thr'])
        
        sequence = []
        for i in range(slength):
            sequence.append(random.choice(ghmm.DNA.listOfCharacters))
        self.tSeq  = ghmm.EmissionSequence(ghmm.DNA, sequence, labelDomain=self.l_domain,labelInput=self.labels)

    def testReadHMMed(self):
        model = ghmm.HMMOpenXML('multexon-4.xml')
        del model
        model = ghmm.HMMOpenXML('test2.xml')
        del model

    def testWriteReadXML(self):
        """
        Test writing from matrices to XML.
        Ignored attributes: tied_to and background.
        """
        self.model.toXML('./discrete.xml')
        model2 = ghmm.HMMOpenXML('./discrete.xml')

        self.label_model.toXML('./model_label.xml')
        model3 = ghmm.HMMOpenXML('./model_label.xml')


# Run ALL tests
if __name__ == '__main__':
    unittest.main()


# Individual test suites for each of the different classes
suiteAlphabet = unittest.makeSuite(AlphabetTests,'test')
suiteEmissionSequence = unittest.makeSuite(EmissionSequenceTests,'test')
suiteSequenceSet = unittest.makeSuite(SequenceSetTests,'test')
suiteDiscreteEmissionHMM = unittest.makeSuite(DiscreteEmissionHMMTests,'test')
suiteBackgroundDistribution = unittest.makeSuite(BackgroundDistributionTests,'test')
suiteStateLabelHMM = unittest.makeSuite(StateLabelHMMTests,'test')
suiteGaussianEmissionHMM = unittest.makeSuite(GaussianEmissionHMMTests,'test')
suiteGaussianMixtureHMM = unittest.makeSuite(GaussianMixtureHMMTests,'test')
suiteXMLIO = unittest.makeSuite(XMLIOTests,'test')

# Call to individual test suites, uncomment to activate as needed.
#runner = unittest.TextTestRunner()
#runner.run(suiteAlphabet)
#runner.run(suiteSequenceSet)
#runner.run(suiteDiscreteEmissionHMM)
#runner.run(suiteBackgroundDistribution)
#runner.run(suiteStateLabelHMM)
#runner.run(suiteGaussianEmissionHMM)
#runner.run(suiteGaussianMixtureHMM)
#runner.run(suiteXMLIO)








