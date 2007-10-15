#!/usr/bin/env python
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

# adjust verbosity level
import logging, sys

log = logging.getLogger("GHMM unit tests")
# creating StreamHandler to stderr
hdlr = logging.StreamHandler(sys.stderr)
# setting message format
fmt = logging.Formatter("%(name)s %(filename)s:%(lineno)d  - %(message)s")
hdlr.setFormatter(fmt)
# adding handler to logger object
log.addHandler(hdlr)
#set unittests log level
log.setLevel(logging.WARN)

#set GHMM log level
ghmm.log.setLevel(logging.WARN)


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
        # remove this assertion because -1 now represents a gap '-'
        # self.assertRaises(KeyError, self.dnaAlphabet.external, -1)
        
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
        #print"\ntestlen ",
        
        self.assertEqual(len(self.binaryAlphabet),2)
        self.assertEqual(len(self.dnaAlphabet),len(self.dna))

        self.assertEqual(len(self.binaryAlphabet),2)
        self.assertEqual(len(self.dnaAlphabet),len(self.dna))
    

class EmissionSequenceTests(unittest.TestCase):
    
    def setUp(self):
        self.i_dom = ghmm.IntegerRange(0,5)
        self.d_dom = ghmm.Float()
        l_domain = ghmm.LabelDomain(['E','R','T'])
        self.i_seq = ghmm.EmissionSequence(self.i_dom,[1,2,0,0,0,3,4])
        self.d_seq = ghmm.EmissionSequence(self.d_dom,[1.3, 2.1, 0.8, 0.1, 0.03, 3.6, 43.3])
        self.labeled = ghmm.EmissionSequence(ghmm.DNA, list('acgttgatgga'),labelDomain=l_domain,
                                            labelInput= ['E','R','T','T','T','E','R','T','T','T','R'])


    def testprint(self):
        #print"\ntestprint ",
        s = "\nEmissionSequence Instance:\nlength 7, weight 1.0:\n1200034"
        self.assertEqual(self.i_seq.verboseStr(),s)
        
        s2 = "\nEmissionSequence Instance:\nlength 7, weight 1.0:\n1.3 2.1 0.8 0.1 0.03 3.6 43.3 "
        self.assertEqual(self.d_seq.verboseStr(),s2)
        
        
    def testattributes(self):
        #print"\ntestattributes ", 
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
    
    def testFileIO(self):
        # print"\ntestFileIO ",
        self.i_seq.write("testdata/es_discrete_testwrite.seq")
        self.d_seq.write("testdata/es_continuous_testwrite.seq")

        discrete_seq   = ghmm.EmissionSequence(self.i_dom, "testdata/es_discrete_testwrite.seq")
        continuous_seq = ghmm.EmissionSequence(self.d_dom, "testdata/es_continuous_testwrite.seq")
        
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
        self.i_seq.setSeqLabel(8)
        l = self.i_seq.getSeqLabel()   
        self.assertEqual(l,8)
        
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
        label = self.labeled.getStateLabel()
        self.assertEqual(label,['E','R','T','T','T','E','R','T','T','T','R'])


class SequenceSetTests(unittest.TestCase):
    
    def setUp(self):
        #print "----------------- Setting up... ---------"
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
        #print"\n----------------- testprint "
        s = "\nNumber of sequences: 5\nSeq 0, length 5, weight 1.0:\n12345\nSeq 1, length 3, weight 1.0:\n030\nSeq 2, length 8, weight 1.0:\n43221111\nSeq 3, length 5, weight 1.0:\n00021\nSeq 4, length 6, weight 1.0:\n111111"
        self.assertEqual(self.i_seq.verboseStr(), s)

        s2 = "\nNumber of sequences: 5\nSeq 0, length 5, weight 1.0:\n1.5 2.3 3.7 4.1 5.1 \nSeq 1, length 3, weight 1.0:\n0.0 3.1 0.7 \nSeq 2, length 8, weight 1.0:\n4.4 3.05 2.0 2.4 1.2 1.8 1.0 1.0 \nSeq 3, length 5, weight 1.0:\n0.4 0.1 0.33 2.7 1.345 \nSeq 4, length 6, weight 1.0:\n1.0 1.0 1.0 1.0 1.0 1.0 "
        self.assertEqual(self.d_seq.verboseStr(), s2)

        # XXX str(self.l_seq)

       
    def testattributes(self):
        #print"\n----------------- testattributes "
        self.assertEqual(len(self.i_seq),5)
        self.assertEqual(self.i_seq.sequenceLength(1),3)
        
        self.assertEqual(len(self.d_seq),5)
        self.assertEqual(self.d_seq.sequenceLength(4),6)
     
    def testgetitem(self):
        #print"\n----------------- testgetitem ",
        s = self.i_seq[2]
        self.assertEqual(len(s),8)
        
        s2 = self.d_seq[4]
        self.assertEqual(len(s2),6)
        
    
    def testweightaccess(self):
        #print"\n----------------- testweightaccess "
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
        #print"\n************************************ testmerge "
        wrong = 4  # wrong argument type to merge
        self.assertRaises(TypeError,self.i_seq.merge,wrong)

        mseq = ghmm.SequenceSet(self.i_alph,[ [1,4,0,4,5,3],[1,2,3,0] ])        
        self.i_seq.merge(mseq)
        self.assertEqual(len(self.i_seq),7)
        s = "\nNumber of sequences: 7\nSeq 0, length 5, weight 1.0:\n12345\nSeq 1, length 3, weight 1.0:\n030\nSeq 2, length 8, weight 1.0:\n43221111\nSeq 3, length 5, weight 1.0:\n00021\nSeq 4, length 6, weight 1.0:\n111111\nSeq 5, length 6, weight 1.0:\n140453\nSeq 6, length 4, weight 1.0:\n1230"
        self.assertEqual(self.i_seq.verboseStr(),s)

        #print self.i_seq

        d_mseq = ghmm.SequenceSet(self.d_alph,[ [7.5,4.0,1.2],[0.4,0.93,3.3,2.54] ])    
        self.d_seq.merge(d_mseq)
        self.assertEqual(len(self.d_seq),7)
        s2 = "\nNumber of sequences: 7\nSeq 0, length 5, weight 1.0:\n1.5 2.3 3.7 4.1 5.1 \nSeq 1, length 3, weight 1.0:\n0.0 3.1 0.7 \nSeq 2, length 8, weight 1.0:\n4.4 3.05 2.0 2.4 1.2 1.8 1.0 1.0 \nSeq 3, length 5, weight 1.0:\n0.4 0.1 0.33 2.7 1.345 \nSeq 4, length 6, weight 1.0:\n1.0 1.0 1.0 1.0 1.0 1.0 \nSeq 5, length 3, weight 1.0:\n7.5 4.0 1.2 \nSeq 6, length 4, weight 1.0:\n0.4 0.93 3.3 2.54 "
        self.assertEqual(self.d_seq.verboseStr(),s2)        
        
    
    def testgetsubset(self):
        #print"\n----------------- testgetsubset "
        i_subseq = self.i_seq.getSubset([2,1,3])
        s = "\nNumber of sequences: 3\nSeq 0, length 8, weight 1.0:\n43221111\nSeq 1, length 3, weight 1.0:\n030\nSeq 2, length 5, weight 1.0:\n00021"
        self.assertEqual(i_subseq.verboseStr(),s)
        self.assertEqual(len(i_subseq),3)
        self.assertEqual(i_subseq.sequenceLength(0),8)
        
        d_subseq = self.d_seq.getSubset([0,4])
        s2 = "\nNumber of sequences: 2\nSeq 0, length 5, weight 1.0:\n1.5 2.3 3.7 4.1 5.1 \nSeq 1, length 6, weight 1.0:\n1.0 1.0 1.0 1.0 1.0 1.0 "
        self.assertEqual(d_subseq.verboseStr(),s2)
        self.assertEqual(len(d_subseq),2)
        self.assertEqual(d_subseq.sequenceLength(0),5)
        
    def testwrite(self):
        #print"\n-----------------testwrite "
        self.i_seq.write("testdata/ghmmunittests_testwrite.seq") 
        self.d_seq.write("testdata/ghmmunittests_testwrite.seq") 

       
    def testlabelaccess(self):
        #print "\n ----------------- testlabelaccess"
        self.i_seq.getSeqLabel(2)   
        l = self.d_seq.getSeqLabel(3)
        self.assertEqual(l,-1)
        self.d_seq.setSeqLabel(3,8)
        l = self.d_seq.getSeqLabel(3)
        self.assertEqual(l,8)

    def testfilereading(self):
        #print "\n ----------------- testfilereading"
        dom = ghmm.IntegerRange(0,12)
        seqs = ghmm.SequenceSetOpen(dom, 'testdata/d_seq.sqd')
        seqs = ghmm.SequenceSetOpen(self.d_alph, 'testdata/test10.sqd')
        seqs = ghmm.SequenceSetOpen(ghmm.Float(), 'testdata/tiny.txt.sqd')

class HMMBaseClassTests(unittest.TestCase):
    def setUp(self):
        A   = [[0.3,0.3,0.4],[0.6,0.1,0.3],[1.0,0.0,0.0]]
        B   = [[0.0,0.5,0.5,0.0],[0.1,0.0,0.8,0.1], [0.0,0.0,0.0,0.0]]
        pi  = [1.0,0.0,0.0]
        tmp = ghmm.HMMFromMatrices(ghmm.DNA, ghmm.DiscreteDistribution(ghmm.DNA), A, B, pi)
        self.model = ghmm.HMM(ghmm.DNA, ghmm.DiscreteDistribution(ghmm.DNA), tmp.cmodel)

    def testpathPosteriorExeption(self):
        self.assertRaises(NotImplementedError, self.model.pathPosterior, [1,2], 34)

    def teststatePosteriorExeption(self):
        self.assertRaises(NotImplementedError, self.model.statePosterior, "sequence", "state", "time")

    def testposteriorExeption(self):
        self.assertRaises(NotImplementedError, self.model.posterior, "sequence")

    def testlogprobExeption(self):
        self.assertRaises(NotImplementedError, self.model.logprob, "emissionSequence", "stateSequence")

    def testbaumWelchExeption(self):
        self.assertRaises(NotImplementedError, self.model.baumWelch, "trainingSequences", "nrSteps", "loglikelihoodCutoff")

    def testbaumWelchSetupExeption(self):
        self.assertRaises(NotImplementedError, self.model.baumWelchSetup, "trainingSequences", "nrSteps")

    def testbaumWelchStepExeption(self):
        self.assertRaises(NotImplementedError, self.model.baumWelchStep, "nrSteps", "loglikelihoodCutoff")
    
    def testbaumWelchDeleteExeption(self):
        self.assertRaises(NotImplementedError, self.model.baumWelchDelete)

    def teststateExeption(self):
        self.assertRaises(NotImplementedError, self.model.state, "stateLabel")

    def testsetEmissionExeption(self):
        self.assertRaises(NotImplementedError, self.model.setEmission, 1, "blah")

    def testtoMatricesExeption(self):
        self.assertRaises(NotImplementedError, self.model.toMatrices)

    def testnormalizeExeption(self):
        self.assertRaises(NotImplementedError, self.model.normalize)

    def testrandomizeException(self):
        self.assertRaises(NotImplementedError, self.model.randomize, "noiseLevel")


class DiscreteEmissionHMMTests(unittest.TestCase):
    def setUp(self):
        log.debug("DiscreteEmissionHMMTests.setUp() -- begin")
        self.A = [[0.3,0.3,0.4],[0.6,0.1,0.3],[1.0,0.0,0.0]]
        self.B = [[0.0,0.5,0.5,0.0],[0.1,0.0,0.8,0.1], [0.0,0.0,0.0,0.0]]
        self.pi = [1.0,0.0,0.0]
        self.model = ghmm.HMMFromMatrices(ghmm.DNA,ghmm.DiscreteDistribution(ghmm.DNA), self.A, self.B, self.pi)
        log.debug("DiscreteEmissionHMMTests.setUp() -- end")

    def test__str__(self):
        log.debug("test__str__ -- begin")
        # we aren't interested in the output but the function should run fine
        str(self.model)

    def testAccessFunctions(self):
        log.debug("testAccessFunctions -- begin")

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
        self.assertEqual(self.model.cmodel.model_type,260)
        self.assertEqual(self.model.getSilentFlag(1),1)
        
        # removing silent state
        self.model.setEmission(1,[0.2,0.2,0.2,0.4])
        emission = self.model.getEmission(1)
        self.assertEqual(emission,[0.2,0.2,0.2,0.4] )  
        self.assertEqual(self.model.cmodel.model_type,260)
        self.assertEqual(self.model.getSilentFlag(1),0)
        
        # removing last silent state
        self.model.setEmission(2,[0.25,0.25,0.25,0.25])
        emission = self.model.getEmission(2)
        self.assertEqual(emission,[0.25,0.25,0.25,0.25])  
        self.assertEqual(self.model.cmodel.model_type,256)
        self.assertEqual(self.model.getSilentFlag(2),0)
    
        # inserting silent state
        self.model.setEmission(2,[0.0,0.0,0.0,0.0])
        emission = self.model.getEmission(2)
        self.assertEqual(emission,[0.0,0.0,0.0,0.0])  
        self.assertEqual(self.model.cmodel.model_type,260)
        self.assertEqual(self.model.getSilentFlag(2),1)

    def testNewXML(self):
        log.debug("testNewXML -- begin")
        model = ghmm.HMMOpenXML('../doc/xml_example.xml')
    
    def getModel(self):
        log.debug("getModel -- begin")
        A  = [[0.3, 0.6,0.1],[0.0, 0.5, 0.5],[0.0,0.0,1.0]]
        B  = [[0.5, 0.5],[0.5,0.5],[1.0,0.0]]
        pi = [1.0, 0.0, 0.0]
        return ghmm.HMMFromMatrices(ghmm.IntegerRange(0,2),
                                    ghmm.DiscreteDistribution(ghmm.IntegerRange(0,2)),
                                    A, B, pi)
        
    def testDel(self):
        """  test for explicit construction and destruction """
        log.debug("testDel -- begin")
        del self.model
        print "===== testing for construction/destruction ===="
        for i in range(100):
            mo = self.getModel()
    
    def testAsMatrices(self):
        log.debug("testAsMatrices -- begin")
        tA,tB,tpi = self.model.toMatrices()
        
        self.assertEqual(self.A,tA)
        self.assertEqual(self.B,tB)
        self.assertEqual(self.pi,tpi)        
            
    def testSample(self):
        log.debug("testSample -- begin")
        seq = self.model.sampleSingle(100, seed=3586662)
        #print seq
        seq2 = self.model.sample(10, 100, seed=3586662)


    def testBaumWelch(self):
        log.debug("testBaumWelch -- begin")
        seq = self.model.sample(100,100,seed=3586662)
        self.assertRaises(NotImplementedError, self.model.baumWelch, seq, 5, 0.01)
        
        self.model.setEmission(2,[0.25,0.25,0.25,0.25])
        self.assertEqual(self.model.cmodel.model_type & 4, 0)
        self.model.baumWelch(seq,5,0.01)
        self.model.baumWelch(seq)
        
    def testViterbi(self):
        log.debug("testViterbi -- begin")
        # Caution with an even number of consecutives 'g'
        # the model can produce two equal probable paths
        
        f = lambda x: round(x, 13)
        g = lambda x: map(f, x)

        seq = ghmm.EmissionSequence(ghmm.DNA, [ 'c','c','g','c','c','g','g','g','g','g','c','g','g','g','c' ])
                                               #[0,  2,  0,  1,  0,  2,  0,  1,  0,  1,  0,  1,  0,  1,  0]
                                               #[0,  2,  0,  1,  0,  2,  0,  1,  0,  1,  0,  1,  0,  1,  0, 1, 0]
                
        result = self.model.viterbi(seq)
        trueResult =  ([0, 2, 0, 1, 0, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0], -19.698557965224637)
        
        self.assertEqual (result[0], trueResult[0])
        self.assertEqual (f(result[1]), f(trueResult[1]))
    
        
        seq2 = ghmm.SequenceSet(ghmm.DNA, [['c','c','g','c','c','g','g','g','g','g','c','g','g','g','c'],
                                           ['c','g','g','g','g','g','c','t','g','c','g','g','t','c','c'],
                                           ['g','g','c','g','c','c','g','c','c','c','c','g','g','g','t'],
                                           ['g','c','c','g','c','c','g','g','g','c','c','c','g','g','g'],
                                           ['g','g','g','g','c','g','c','a','g','g','c','c','g','g','g'],
                                           ['g','g','c','g','g','t','g','c','c','c','c','g','g','a','a'],
                                           ['c','g','g','g','c','c','t','g','c','g','c','g','a','g','g'],
                                           ['c','c','g','g','g','g','g','g','g','c','g','c','g','c','g'],
                                           ['g','g','g','g','c','c','g','g','g','c','g','c','g','g','g'],
                                           ['g','c','t','c','g','g','a','g','g','c','a','g','g','g','g']])
        
        
        path2 = self.model.viterbi(seq2)
        truePath2 = ([[0, 2, 0, 1, 0, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
                      [0, 1, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 0],
                      [0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 2, 0, 1, 0, 1],
                      [0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 2, 0, 1, 0, 1],
                      [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1],
                      [0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2, 0, 2, 0, 1, 0, 1, 1],
                      [0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 1],
                      [0, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                      [0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                      [0, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]],
                      [-19.698557965224637, -23.857441048584306, -22.77125127991475,
                       -20.691809738234902, -21.595677950110499, -27.846425095148575,
                       -24.262906156692473, -19.516236408430682, -19.5162364084307,
                       -25.754561033470189])
                
        self.assertEqual (path2[0], truePath2[0])    
        self.assertEqual (g(path2[1]), g(truePath2[1]))    
        

    def testLoglikelihood(self):
        log.debug("testLoglikelihood -- begin")
        seq = self.model.sampleSingle(100,seed=3586662)
        logp = self.model.loglikelihood(seq)
        self.assert_(logp-93.1053904716 < 10^-8, "Different results in loglikelihood ")

    def testLogProb(self):
        log.debug("testLogProb -- begin")
        seq = self.model.sampleSingle(15,seed=3586662)
        path,logp = self.model.viterbi(seq)
        logp = self.model.logprob(seq,path)
        self.assert_(logp - 22.4303246929 < 10^-8, "Different results in logprob ")
        
    def testFoBa(self):
        log.debug("testFoBa -- begin")
        #print self.model
        seq = ghmm.EmissionSequence(self.model.emissionDomain,
                                    ['g','g','g','c','t','g','g','c','g','g',
                                     'g','c','g','g','g','c','c','c','g','c',
                                     'g','g','c','c','g','c','g','c','c','c',
                                     'g','c','g','c','g','g','c','t','c','c'])

        (alpha,scale) = self.model.forward(seq)
       
        
        f = lambda x: round(x, 13)
        g = lambda x: map(f, x)
        
        talpha = [[0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.7142857142857143, 0.0, 0.2857142857142857], [0.0, 0.76923076923076916, 0.23076923076923075], [0.61307901907356943, 0.10899182561307905, 0.27792915531335155], [0.46113149992850677, 0.27262761546160807, 0.26624088460988515], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.48775986167826385, 0.24395091819263898, 0.26828922012909723], [0.71428571428571419, 0.0, 0.2857142857142857], [0.43640897755610969, 0.29925187032418948, 0.26433915211970072], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.48775986167826385, 0.24395091819263898, 0.26828922012909723], [0.71428571428571419, 0.0, 0.2857142857142857], [0.7142857142857143, 0.0, 0.2857142857142857], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.7142857142857143, 0.0, 0.2857142857142857], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.7142857142857143, 0.0, 0.28571428571428575], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.7142857142857143, 0.0, 0.28571428571428575], [0.43640897755610975, 0.29925187032418954, 0.26433915211970077], [0.50453092851201675, 0.22588976929475116, 0.26957930219323206], [0.7142857142857143, 0.0, 0.2857142857142857], [0.0, 0.76923076923076916, 0.23076923076923075], [0.71428571428571419, 0.0, 0.2857142857142857], [0.7142857142857143, 0.0, 0.2857142857142857]]
        self.assertEqual(map(g, alpha),map(g, talpha))
        
        tscale =   [0.69999999999999996, 0.57285714285714284, 0.56965087281795512, 0.38953070962658143, 0.027857142857142858, 0.56461538461538452, 0.57168937329700276, 0.39770983270578136, 0.57285714285714284, 0.56965087281795512, 0.57043689532898478, 0.39269141068371188, 0.57285714285714284, 0.56965087281795501, 0.57043689532898478, 0.39269141068371188, 0.34999999999999998, 0.34999999999999998, 0.57285714285714284, 0.40236907730673316, 0.57285714285714284, 0.56965087281795512, 0.38953070962658143, 0.34999999999999998, 0.57285714285714284, 0.40236907730673316, 0.57285714285714284, 0.40236907730673316, 0.34999999999999998, 0.34999999999999998, 0.57285714285714284, 0.40236907730673316, 0.57285714285714284, 0.40236907730673316, 0.57285714285714284, 0.56965087281795512, 0.38953070962658143, 0.027857142857142858, 0.48461538461538456, 0.34999999999999998]
        self.assertEqual(map(f, scale), map(f,tscale))
        
        beta = self.model.backward(seq,scale)
        tbeta = [[0.99999999999999944, 0.93777288282264037, 0.90665932423396078], [1.0387725400509094, 0.87202814099718418, 0.78865594147032159], [0.89851709082326969, 1.1552362596299182, 1.2835958440332425], [0.99999999999999956, 0.3333333333333332, 0.0], [0.99018797150167415, 0.92857142857142794, 0.89776315710630483], [1.0137817804861966, 0.8510489133365684, 0.76968247976175441], [0.88003858898536058, 1.1314781858383207, 1.2571979842648009], [0.999999999999999, 0.91297745742272951, 0.86946618613409465], [0.99615983039934863, 0.93417167590571015, 0.90317759865889091], [1.028991814771324, 0.86381742368004855, 0.7812302281344109], [0.89128509174829584, 1.1459379751049517, 1.2732644167832796], [0.99999999999999944, 0.91297745742272962, 0.86946618613409488], [0.99615983039934886, 0.93417167590571037, 0.90317759865889113], [1.0289918147713242, 0.86381742368004866, 0.78123022813441101], [0.89128509174829595, 1.1459379751049521, 1.2732644167832801], [0.99999999999999967, 1.2857142857142854, 1.4285714285714282], [1.0, 1.2857142857142856, 1.4285714285714286], [1.0, 0.83947939262472882, 0.75921908893709322], [0.86984815618221256, 1.1183762008057019, 1.2426402231174465], [1.0, 0.93777288282264093, 0.90665932423396134], [1.03877254005091, 0.87202814099718462, 0.78865594147032214], [0.89851709082327025, 1.1552362596299188, 1.2835958440332431], [1.0000000000000002, 1.285714285714286, 1.428571428571429], [1.0000000000000002, 0.83947939262472893, 0.75921908893709344], [0.86984815618221256, 1.1183762008057019, 1.2426402231174465], [1.0, 0.83947939262472882, 0.75921908893709322], [0.86984815618221256, 1.1183762008057019, 1.2426402231174465], [1.0, 1.2857142857142856, 1.4285714285714286], [1.0, 1.2857142857142856, 1.4285714285714286], [1.0, 0.83947939262472882, 0.75921908893709322], [0.86984815618221256, 1.1183762008057019, 1.2426402231174465], [1.0, 0.83947939262472882, 0.75921908893709322], [0.86984815618221256, 1.1183762008057019, 1.2426402231174465], [1.0, 0.93777288282264093, 0.90665932423396134], [1.03877254005091, 0.87202814099718462, 0.78865594147032214], [0.89851709082327025, 1.1552362596299188, 1.2835958440332431], [1.0000000000000002, 0.33333333333333343, 0.0], [0.72222222222222221, 0.92857142857142849, 1.0317460317460319], [1.0, 1.2857142857142856, 1.4285714285714286], [1.0, 1.0, 1.0]]
        
        self.assertEqual (map(g, beta),map(g, tbeta))
                
        #testing forward and backward log probabilities
        self.assertEqual (f(self.model.loglikelihood(seq)),
                          f(self.model.backwardTermination (seq, beta, scale)))

    def testTiedStates(self):
        log.debug( "testTiedStates -- begin")
        f = lambda x: round(x,15)
        t = (-1,1,1)
        self.model.setTieGroups(t)

        self.model.updateTieGroups()
        em2 = map(f,self.model.getEmission(2))
        self.assertEqual(em2, [0.0, 0.0, 0.0, 0.0])

        self.model.setEmission(2,[0.2,0.2,0.2,0.4])
        self.model.updateTieGroups()
        em0 = map(f,self.model.getEmission(0))
        self.assertEqual(em0, [0.0,0.5,0.5,0.0])
        em2 = map(f,self.model.getEmission(2))
        self.assertEqual(em2, [0.15, 0.1, 0.5, 0.25])
        log.debug("testTiedStates -- end")

    def testNormalization(self):
        log.debug("testNormalization")
        self.model.setInitial(0, 2.0)
        self.model.normalize()
        self.assertEqual(1.0, self.model.getInitial(0))
        
class BackgroundDistributionTests(unittest.TestCase):
    " Tests for background distributions "

    def setUp(self):
        self.sigma = ghmm.Alphabet(['rot','blau','gruen','gelb'])

        self.model = ghmm.HMMFromMatrices(self.sigma,ghmm.DiscreteDistribution(self.sigma),
                       [[0.3,0.3,0.4],[0.6,0.1,0.3],[1.0,0.0,0.0]],
                       [[0.0,0.5,0.5,0.0],[0.1,0.0,0.8,0.1],
                       [0.25,0.25,0.25,0.25, 0.0,0.5,0.5,0.0, 0.1,0.0,0.8,0.1, 0.1,0.35,0.3,0.25 ]],
                       [1.0,0,0])

        self.bg = ghmm.BackgroundDistribution(self.sigma, [[0.2,0.3,0.1,0.4],
                       [0.1,0.2,0.4,0.3, 0.2,0.3,0.1,0.4, 0.25,0.25,0.25,0.25, 0.0,0.5,0.5,0.0 ]]
                                              )

    def test__str__(self):
        # we aren't interested in the output but the function should run fine
        str(self.model)
        
    def testprint(self):
        #print "*** testprint"
        s = self.bg.verboseStr()
        ts = "BackgroundDistribution instance:\nNumber of distributions: 2\n\n<Alphabet:['rot', 'blau', 'gruen', 'gelb']>\nDistributions:\n  Order: 0\n  1: [0.20000000000000001, 0.29999999999999999, 0.10000000000000001, 0.40000000000000002]\n  Order: 1\n  2: [0.10000000000000001, 0.20000000000000001, 0.40000000000000002, 0.29999999999999999]\n"
        self.assertEqual(s,ts)

    def testmodelbackgroundaccessfunctions(self):
        
        #print "***  testmodelbackgroundaccessfunctions"
        
        self.model.setBackground(self.bg, [0,-1,1])
        # deleting background
        del(self.bg)
        s = self.model.background.verboseStr()
        ts = "BackgroundDistribution instance:\nNumber of distributions: 2\n\n<Alphabet:['rot', 'blau', 'gruen', 'gelb']>\nDistributions:\n  Order: 0\n  1: [0.20000000000000001, 0.29999999999999999, 0.10000000000000001, 0.40000000000000002]\n  Order: 1\n  2: [0.10000000000000001, 0.20000000000000001, 0.40000000000000002, 0.29999999999999999]\n"
        self.assertEqual(s,ts)

    def testapplybackground(self):
        self.model.setBackground(self.bg,[0, -1, 1])
        self.model.applyBackground([0.1, 0.2, .3])
        #print self.model
        
        f = lambda x: round(x,15)
        e1 = map(f, self.model.getEmission(0)) 
        e2 = map(f, self.model.getEmission(1))
        e3 = map(f, self.model.getEmission(2))
                
        self.assertEqual(e1, [0.02, 0.48, 0.46, 0.04])
        self.assertEqual(e2, [0.1,  0.0,  0.8,  0.1])
        self.assertEqual(e3, [0.205, 0.235, 0.295, 0.265, 0.06, 0.44, 0.38, 0.12, 0.145, 0.075, 0.635, 0.145, 0.07, 0.395, 0.36, 0.175])
        

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

    def test__str__(self):
        # we aren't interested in the output but the function should run fine
        str(self.model)

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
        self.assertEqual(ghmmwrapper.int_array_getitem(self.model.cmodel.silent,1),1)

        # removing silent state
        self.model.setEmission(1,[0.2,0.2,0.2,0.4])
        emission = self.model.getEmission(1)
        self.assertEqual(emission,[0.2,0.2,0.2,0.4] )  
        #print "model_type = ",self.model.cmodel.model_type
        self.assertEqual(self.model.cmodel.model_type & 4,0)
        self.assertEqual(self.model.isSilent(1), False)

        # inserting silent state
        self.model.setEmission(0,[0.0,0.0,0.0,0.0])
        emission = self.model.getEmission(0)
        self.assertEqual(emission,[0.0,0.0,0.0,0.0])  
        self.assertEqual(self.model.cmodel.model_type & 4,4)
        self.assertEqual(ghmmwrapper.int_array_getitem(self.model.cmodel.silent,0),1)

        # label access
        labels = self.model.getLabels()
        self.assertEqual(labels,['fst','scd','thr'])
        self.model.setLabels(['fst','thr','fst'])
        labels = self.model.getLabels()
        self.assertEqual(labels, ['fst','thr','fst']) 

    def testonelabelcomparebackward(self):
        model = self.oneModel(['One']*11)

        # backward and labeled backward use numerical different algorithms
        # to be changed
        f = lambda x: round (x, 12)
        g = lambda x: map (f, x)
        
        labelSequence      = self.labels
        (alpha, scale)     = model.forward( self.tSeq)
        (b_beta)           = model.backward( self.tSeq, scale)
        (bl_logp, bl_beta) = model.backwardLabels( self.tSeq, labelSequence, scale)

        #compare beta matrizes from backward and backwardLabels (all states share one label)
        self.assertEqual (map(g, b_beta), map(g, bl_beta))


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
                if model2.labelDomain.internal(labelSequence[i]) == ghmmwrapper.int_array_getitem(model2.cmodel.label, j):
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
        model2  = self.oneModel(self.allLabels)
        
        labelSequence = self.allLabels*4

        sequence = []
        for i in range(len(labelSequence)):
            sequence.append(random.choice(ghmm.DNA.listOfCharacters))

        Seq = ghmm.EmissionSequence(ghmm.DNA, sequence, self.l_domain, labelSequence)

        (logp, alpha, scale) =  model2.forwardLabels(Seq, labelSequence)

        #check if the beta matrix is 0 or different from 0 at the appropriate entries 
        for i in range(len(alpha)):
            i = len(alpha)-i-1
            for j in range(len(alpha[i])):
                if model2.labelDomain.internal(labelSequence[i]) == ghmmwrapper.int_array_getitem(model2.cmodel.label, j):
                    self.assertNotEqual(alpha[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(alpha[i][j]) )
                else:
                    self.assertEqual(alpha[i][j], 0.0, "Zeichen: " + str(i) + ", State: " + str(j)
                                        + ", value: " + str(alpha[i][j]))


    def testkbest(self): 
        seq = ghmm.EmissionSequence(self.model.emissionDomain,
                                    ['a','c','g','t','t','a','a','a','c','g',
                                     't','g','a','c','g','c','a','t','t','t'],
                                    self.model.labelDomain,
                                    ['fst', 'scd', 'thr', 'thr', 'thr', 'thr', 'scd',
                                     'scd', 'thr', 'thr', 'scd', 'thr', 'scd', 'thr',
                                     'thr', 'thr', 'scd', 'fst', 'scd', 'fst'])


        path = self.model.kbest(seq)
        self.assertEqual(path,(['fst', 'thr', 'thr', 'scd', 'fst', 'scd', 'fst', 'scd',
                                'thr', 'thr', 'fst', 'thr', 'thr', 'thr', 'thr', 'thr',
                                'scd', 'fst', 'scd', 'fst'], -35.735009627142446)) 


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
        self.model.baumWelchLabels(seq,5,0.01)   
    
   
    def testviterbilabels(self):
        seq = ghmm.EmissionSequence(ghmm.DNA, ['a','c','g','t','t','a','a','a','c','g','t','g','a','c','g','c','a','t','t','t'])
        
        p = self.model.viterbiLabels(seq)

        self.assertEqual(p[0], ['fst', 'thr', 'thr', 'scd', 'fst', 'scd', 'fst', 'scd', 'thr', 'thr', 'fst', 'thr', 'thr', 'thr', 'thr', 'thr', 'scd', 'fst', 'scd', 'fst'])
        self.assertEqual(round(p[1],14) ,round(-39.893892710502115,14))


    # TO DO: testing XML-file read
   

class GaussianEmissionHMMTests(unittest.TestCase):

    def setUp(self):
        # print"setUp"
        F = ghmm.Float()
        self.A = [[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]]
        self.B = [[0.0,1.0],[-1.0,0.5], [1.0,0.2]]
        self.pi = [1.0,0.0,0.0]
        self.model = ghmm.HMMFromMatrices(F,ghmm.GaussianDistribution(F), self.A, self.B, self.pi)

    def test__str__(self):
        # we aren't interested in the output but the function should run fine
        str(self.model)

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
        #print "** GaussianMixtureHMMTests **"

    def test__str__(self):
        # we aren't interested in the output but the function should run fine
        str(self.model)

    def testSMO(self):
        model = ghmm.HMMOpenSMO('testdata/tiny.smo')

    def testNewXML(self):
        model = ghmm.HMMOpenXML('../doc/xml_cont_example.xml')

    def testMultipleTransitionClasses(self):
        model = ghmm.HMMOpenXML('testdata/xml_cont_multiple.xml')
        state = model.cmodel.getState(0)
        self.assertEqual(state.getOutProb(0, 0), state.getOutProb(0))
        self.assertEqual(state.getOutProb(0, 0), 0.1)
        self.assertEqual(state.getOutProb(0, 1), 0.2)

    def testcomponentfixing(self):
        #print "testcomponentfixing"
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
        #print"\ntesttomatrices "
        tA,tB,tpi = self.model.toMatrices()
        
        self.assertEqual(self.A,tA)
        self.assertEqual(self.B,tB)
        self.assertEqual(self.pi,tpi)    
        
        
    def testsample(self):
        #print"\ntestsample "
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
        #print "** XMLIOTests **"
        #print self.tSeq

    def testReadHMMed(self):
        model = ghmm.HMMOpenXML('testdata/multexon-4.xml')
        del model
        model = ghmm.HMMOpenXML('testdata/test2.xml')
        del model

    def testWriteReadXML(self):
        """
        Test writing from matrices to XML.
        Ignored attributes: tied_to and background.
        """
        #self.model.toXML('testdata/discrete.xml')
        model2 = ghmm.HMMOpenXML('testdata/discrete.xml')

        #self.label_model.toXML('testdata/model_label.xml')
        model3 = ghmm.HMMOpenXML('testdata/model_label.xml')

########### PAIR HMM TESTS ##############

class ComplexEmissionSequenceTests(unittest.TestCase):
    
    def setUp(self):
        i_alph = ghmm.IntegerRange(0,5)
        d_alph = ghmm.Float()
        self.seq = ghmm.ComplexEmissionSequence([i_alph, ghmm.DNA, d_alph],
                                                [[1,2,0,0,0,3,4],
                                                 ['a','t','g','c','t','g','c'],
                                                 [1.3, 2.1, 0.8, 0.1, 0.03, 3.6, 43.3]])
        
    def testprint(self):
        # print"\ntestprint ",
        s = ("ComplexEmissionSequence (len=7, discrete=2, continuous=1)\n" +
             "1200034\n" +
             "atgctgc\n" +
             "1.3,2.1,0.8,0.1,0.03,3.6,43.3\n")

        self.assertEqual(self.seq.verboseStr(),s)
        
        
    def testattributes(self):
        # print"\ntestattributes ", 
        self.assertEqual(self.seq.cseq.number_of_alphabets,2)
        self.assertEqual(self.seq.cseq.number_of_d_seqs,1)    
        self.assertEqual(self.seq.cseq.length,7)    
        self.assertEqual(len(self.seq),7)    
    
    def testitemaccess(self):
        # print"\ntestitemaccess ",
        b = self.seq.getInternalDiscreteSequence(0)
        self.assertEqual(b[5], 3)    
        
        b2 = self.seq.getInternalContinuousSequence(0)
        self.assertEqual(b2[1],2.1)    

    def testerrors(self):
        pass


# Run ALL tests (comment out to deactivate)
#if __name__ == '__main__':
#    unittest.main()


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
runner = unittest.TextTestRunner()
#runner.run(suiteAlphabet)
#runner.run(suiteEmissionSequence)
#runner.run(suiteSequenceSet)
#runner.run(suiteDiscreteEmissionHMM)
#runner.run(suiteBackgroundDistribution)
#runner.run(suiteStateLabelHMM)
#runner.run(suiteGaussianEmissionHMM)
#runner.run(suiteGaussianMixtureHMM)
#runner.run(suiteXMLIO)

