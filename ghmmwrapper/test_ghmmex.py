from new_ghmmex import *
#from gendata import *
from GraphUtil import SaveCATBoxGraph, OpenCATBoxGraph


def get_hmmer_raw_string(fileName):
	"returns the first model in the flat file as a raw string"
	f = open(fileName,"r")

	string = ""
	res = re.compile("^//")
	
	print "reading first model in " + str(fileName) + "."
	for line in f.readlines():
		string = string + line
		match = res.match(line)
		if match:
			break
	f.close()
	return string


gsl_rng_init() # Init random number generator
time_seed()

#string = get_hmmer_raw_string("simple2.hmm")
lisalph = ['A','C','G','T']

#[matrans,alphabet,maems,pi,acc] = hmmer2ghmmex(hmmer(string),return_ghmm=0)

print "*** TEST:  model.c ***\n"

matrans = [[0.0,1.0,0],[0.5,0.0,0.5],[1.0,0.0,0.0]]
lisalph = Nucleotide # ["A","C","G","T"]
maems = [[1.0,0.0,0.0,0.0],[0.0,0.5,0.5,0.0],[0.0,0.0,1.0,0.0]]
pi = [1.0,0,0]

print "generating model..."
m = hmm_model(matrans,lisalph,maems,pi,"test_model")
#print m
#m.test()
print "generating sequence..."
s = m.generate_sequence(0,10,1)
print s

test = s.onesequence(0)
#print test
print "\nviterbi algorithm...",
(s,p) = m.make_viterbi(test)
print "\n" + str(s),p

print "\nforward algorithm...",
z = m.make_forward(test)
print "\nForward result: " + str(z) 

####################################################################################

print "\n*** TEST:  sdmodel.c ***\n"

sdmatrans = [[[0.0,1.0,0,0],[0.5,0.0,0.5,0],[1.0,0.0,0.0,0],[0,0,0,1] ],[[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1]] ]
sdlisalph = Nucleotide # ["A","C","G","T"]
sdmaems = [[1.0,0.0,0.0,0.0],[0.0,0.5,0.5,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0]]
sdpi = [1.0,0,0,0]

#matrans = [[[0.0,1.0,0],[0.5,0.0,0.5],[1.0,0.0,0.0]]]
#lisalph = Nucleotide # ["A","C","G","T"]
#maems = [[1.0,0.0,0.0,0.0],[0.0,0.5,0.5,0.0],[0.0,0.0,1.0,0.0]]
#pi = [1.0,0,0]

#m2 = hmm_sdmodel(matrans,lisalph,maems,pi)
print "generating sdmodel..."
m2 = hmm_sdmodel(sdmatrans, sdlisalph, sdmaems, sdpi)
# print m2

print "generating sequence..."
s2 = m2.generate_sequence(0,10,1)
test2 = s2.onesequence(0)
print s2


#print "\nviterbi algorithm...",
#(s,p) = m2.make_viterbi(test2)
#print "\n" + str(s),p

#z2 = m2.make_forward(test2)
#print "\nForward result: " + str(z2) 

####################################################################################

print "\n*** TEST:  smodel.c ***\n"
print "reading continous model from file"
sm = SHMM("test2.smo")

print sm

cseq = sm.generate_sequence(0,10,10,6,100)
print cseq



print_sequence_d("seq_test.sqd",cseq.seq_c,0)
