# Python-wrapper for the ghmm library
# Author       : Benjamin Georgi

from ghmmwrapper import *
import math
import os.path
from modhmmer import *
from random import *

#constants
Nucleotide = ['a','c','g','t']
Amino = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
kNotSpecified = 0
kSilentStates = 4

gsl_rng_init() #Init random number generator


def sum_vector(v):
	"sums the components of vector 'v'"
	s = 0
	for i in range(len(v)):
		s=s+v[i]
	return s


def plist2intarray(carr, liste, nstates):
	"converts python list to c int array"
	for i in range(0,nstates):
		set_arrayint(carr,i,liste[i])

def get_col(mat,index):
	"returns <index>th column from matrix <mat>"
	lis = []
	for i in range(len(mat)):
		lis.append(mat[i][index])
	return lis

def get_3d_col(mat,index):
	"returns all columns of index 'index' in the 3d matrix 'mat'"
	cos = len(mat)
	lis = []
	for c in range(cos):
		lis.append([])
		for i in range(len(mat[0])):
			lis[c].append(mat[c][i][index])
	return lis



def list2arrayd(lisems):
	"converts python list to C double array"
	arrems = double_array(len(lisems))
	for i in range(len(lisems)):
		set_arrayd(arrems,i,lisems[i])
	return arrems

def extract_out(lisprobs):
	"initializes C double arrays state-> out_a from a row of the transition matrix"
	lis = []
	for i in range(len(lisprobs)):
		if lisprobs[i]!=0:
			lis.append(i)
	trans_id = int_array(len(lis))
	trans_prob = double_array(len(lis))
	for i in range(len(lis)):
		set_arrayint(trans_id,i,lis[i])
		set_arrayd(trans_prob,i,lisprobs[lis[i]])
	return [len(lis),trans_id,trans_prob]

def extract_out_probs(lisprobs,cos):
	"extract outgoing transitions from state probability list in an sdmodel"

	lis = []
	# parsing indixes belonging to postive probabilites
	for j in range(cos):
		for i in range(len(lisprobs[0])):
			if lisprobs[j][i] != 0 and i not in lis:
				lis.append(i)
	# print "lis: ", lis			


	trans_id   = int_array(len(lis))
	probsarray = double_2d_array(cos, len(lis)) # C-function
	
	# creating list with positive probabilities
	for k in range(cos):
		for j in range(len(lis)):
			set_2d_arrayd(probsarray,k,j, lisprobs[k][lis[j]])

	trans_prob = twodim_double_array(probsarray, cos, len(lis)) # python CLASS
	
	#print trans_prob
	
	# initializing c state index array
	for i in range(len(lis)):
		set_arrayint(trans_id,i,lis[i])
	return [len(lis),trans_id,trans_prob]

def descale(alpha,scale,t,n):
	"wrapps sdfoba_descale"
	print "make alpha"
	alpha_new = double_2darray(t,n)
	for i in range(t):
		for j in range(n):
			set_2d_arrayd(alpha_new,i,j,0.0)
	print "run descale"
	sdfoba_descale(alpha,scale,t,n,alpha_new)
	print "return alpha"
	return alpha_new


class twodim_double_array:
	"Two-dimensional C-Double Array"

	def __init__(self,array, rows, columns, rowlabels=None, columnlabels=None):
		"Constructor"
		self.array = array
		self.rows = rows
		self.columns = columns
		self.size = (rows,columns)
		self.rowlabels =rowlabels
		self.columnlabels = columnlabels

	def __getitem__(self,index):
		"defines twodim_double_array[index[0],index[1]]"
		return get_2d_arrayd(self.array,index[0],index[1])

	def __setitem__(self,index,value):
		"defines twodim_double_array[index[0],index[1]]"
		set_2d_arrayd(self.array,index[0],index[1],value)

	def __str__(self):
		"defines string representation"
		strout = "\n"
		if (self.columnlabels is not None):
			for k in range(len(self.columnlabels)):
				strout+="\t"
				strout+= str(self.columnlabels[k])
			strout += "\n"
		for i in range(self.rows):
			if (self.rowlabels is not None):
				strout += str(self.rowlabels[i])
			strout += "\t"
			for j in range(self.columns):
				strout += "%2.4f" % self[i,j]
				strout += "\t"
			strout += "\n"
		return strout

def remove_state(matrix, index):
	""" Deletes row and column 'index' from matrix 'matrix', this is equivalent to
	removing state 'index' from the HMM defined by 'matrix'."""
	for i in range(len(matrix)):
		del(matrix[i][index])
	del(matrix[index])	
	return matrix

def hmmer2ghmmex(hmmer,return_ghmm=1):
	"converts a MODHMMER object into matrices to initialize GHMMEX object"
	#trans_mat = hmmer.matrans
	emiss_mat = []
	# intitializing pi vector
	pi = []
	pi.append(1) # always starting in B state
	for i in range( (3 * hmmer.n)+1 ): 
		pi.append(0)
	
	if (hmmer.m == 4):
		alphabet = Nucleotide
		silent = [0,0,0,0]
		equal = [0.25,0.25,0.25,0.25]
	if (hmmer.m == 20):
		alphabet = Amino
		silent = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ]
		equal = [1.0/20.0, 1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0,1.0/20.0]
	
	# conversion of the HMMER emission matrices into GHMMEX Format
	# emmission probs in HMMER: [match=0/insert=1][state][emission-letter]
	# order of states in HMMER transition matrix: N B E J C T M1 I1 D1 M2 I2 D2 ... Mn In Dn
	#emiss_mat.append(silent)   #(equal) 				# N state 
	emiss_mat.append(silent) 			    # B state (silent)
	emiss_mat.append(silent)   				# E state (silent)
	emiss_mat.append(equal)     #(silent)  				# J state 
	#emiss_mat.append(silent)   #(equal)  				# C state 
	#emiss_mat.append(silent)   			  	# T (silent)
	for ind1 in range(hmmer.n):
		emiss_mat.append(hmmer.maems[0][ind1])  # M state
		#silent_states.append(0)
		if (ind1 != hmmer.n-1):
			emiss_mat.append(hmmer.maems[1][ind1]) # I state
		emiss_mat.append(silent) #(silent)	 # D state (silent)
		
	# removing unnecessary states 
	tmp = hmmer.matrans
	tmp = remove_state(tmp,0) # removing N state, matrix indices shifted
	tmp = remove_state(tmp,3) # removing C state, matrix indices shifted
	tmp = remove_state(tmp,3) # removing T state
	hmmer.matrans = tmp
	
	hmmer.matrans[2][0] = 1.0  #1 - alg_dep_self # J->B
	hmmer.matrans[2][2] = 0.0  #alg_dep_self # J->J
		
	hmmer.matrans[1][2] = 1 # E->J
	
	if return_ghmm == 1:
		model = hmm_model(hmmer.matrans,alphabet,emiss_mat,pi,hmmer.acc)
		return model
	else:
		return [hmmer.matrans,alphabet,emiss_mat,pi,hmmer.acc]

class intseq:
	"Python wrapper for one C-Integer-Sequence"
	
	def __init__(self,seq,seqlength,seqname=None, alphabet=None):
		"Constructor"	
		self.seq = seq
		self.length = seqlength
		self.name = seqname
		self.alphabet = alphabet

	def __len__(self):
		"Defines lenght of wrapped int array."
		return self.length

	def char(self,x):
		"Returns char belonging to int x as determined by 'self.alphabet'."
		if (self.alphabet is not None): return self.alphabet[x]
		else: return str(x)
		
	def __getitem__(self,index): 
		"Defines 'intseq[index]'."
		assert index >= 0
		assert index < self.length
		return get_arrayint(self.seq,index)
		
	def __str__(self):
		"Defines string representation."
		strout = ""
		if (self.name is not None): strout += str(self.name)+": \n"
		for j in range(self.length): # print each symbol
			strout += self.char(self[j])
		return strout
		
class doubleseq:
	"Python wrapper for one C double arrray."
	
	def __init__(self,seq,seqlength,seqname=None):
		"Constructor"
		self.seq = seq
		self.length = seqlength
		self.name = seqname

	def __len__(self):
		return self.length	
		
	def __getitem__(self,index): 
		"Defines inzseq[index]."
		assert index >= 0
		assert index < self.length
		return get_arrayd(self.seq,index)
	
	def __str__(self):
		"Defines string representation."	 
		strout = ""
		if (self.name is not None): strout += str(self.name)+": \n"
		for j in range(self.length): # print each symbol
			strout += str(get_arrayd(self.seq,j)) + " "
		return strout	
		
		
class smo_array:
	""" Wrapper for a C array of smodel structs."""
	def __init__(self,length):
		"Constructor"
		self.smo_array = smodel_array(length)

	def __getitem__(self,index):
		"Defines 'smo_array[index]'."
		assert index >= 0
		return get_smodel_ptr(self.smo_array, index)
		
	def __setitem__(self,index,ptr):
		"Assigns 'smo_array[index]' = 'ptr'."
		set_smodel_ptr(self.smo_array,ptr, index)
	

class sequence:
	"Python wrapper for ghmm sequence_t struct"

	def __init__(self,seq_c,alphabet=None):
		"Constructor."	
		self.seq_c = seq_c
		self.alphabet = alphabet
	def char(self,x):
		if (self.alphabet is not None): return self.alphabet[x]
		else: return str(x)
	 
	def __getitem__(self,index1,index2=None): 
		"""Defines 'sequence[index1][index2]'. For 'sequence[index1]' the whole sequence 
		   'index1' is returned as an intseq object initialized with a pointer to the sequence
		   in 'sequence' itself.
		"""
		if index2 == None:
			return intseq(get_row_pointer_int(self.seq_c.seq,index1),get_arrayint(self.seq_c.seq_len,index1),"Seq. "+ str(index1))
		else:
			return get_2d_arrayint(self.seq_c.seq,index1,index2)
	 
	def __str__(self):
		"Defines string representation."
		seq = self.seq_c
		strout =  "Number of sequences: " + str(seq.seq_number)
		for i in range(seq.seq_number):
			strout += "\nSeq " + str(i)+ ", length " + str(get_arrayint(seq.seq_len,i))+ ", weight " + str(get_arrayd(seq.seq_w,i))  + ":\n"
			for j in range(get_arrayint(seq.seq_len,i) ):
				strout += self.char(self[i][j])
		return strout

	#def onesequence(self,seqnumber):
		#"Extracts sequence seqnumber from Sequences-Object and returns a intseq object."
		#seq = self.seq_c
		#myseqlength = get_arrayint(seq.seq_len,seqnumber)
		#myseq = int_array(myseqlength)
		#for i in range(myseqlength):
		#	set_arrayint(myseq,i,get_2d_arrayint(seq.seq,seqnumber,i))
		#name = "Sequence "+str(seqnumber)
		#return intseq(myseq,myseqlength,name, self.alphabet)
		
class sequence_d:
	"Python wrapper for ghmm sequence_d_t struct"

	def __init__(self,seq_c):
		"Constructor."	
		self.seq_c = seq_c
	 
	def __getitem__(self,index1,index2=None): 
		"""Defines 'sequence[index1][index2]'. For 'sequence[index1]' the whole sequence 
		   'index1' is returned as a doubleseq object initialized with a pointer to the sequence
		   in 'sequence' itself.
		"""
		if index2 == None:
			return doubleseq(get_row_pointer_d(self.seq_c.seq,index1),get_arrayint(self.seq_c.seq_len,index1),"Seq. "+ str(index1))
		else:
			return get_2d_arrayd(self.seq_c.seq,index1,index2)
	 
	def __str__(self):
		"Defines string representation."
		seq = self.seq_c
		strout =  "\nNumber of sequences: " + str(seq.seq_number)
		for i in range(seq.seq_number):
			strout += "\nSeq " + str(i)+ ", lenght " + str(get_arrayint(seq.seq_len,i))+ ":\n"
			for j in range(get_arrayint(seq.seq_len,i) ):
				strout += str( self[i][j] ) + " "
		return strout		
	
	#def onesequence(self,seqnumber):
		#"Extracts sequence seqnumber from sequence_d object and returns a doubleseq object."
		#seq = self.seq_c
		#myseqlength = get_arrayint(seq.seq_len,seqnumber)
		#myseq = double_array(myseqlength)
		#for i in range(myseqlength):
		#	set_arrayd(myseq,i,get_2d_arrayd(seq.seq,seqnumber,i))
		#name = "Sequence "+str(seqnumber)
		#return doubleseq(myseq,myseqlength,name)		
		

def dummy_2args(a,b):
	"dummy function, 2 arguments"
	print "  dummy 2"
	pass

def dummy_4args(a,b,c,d):
	"dummy function, 4 arguments"
	print "  dummy 4"
	pass

def dummy_5args(a,b,c,d,e):
	"dummy function, 5 arguments"
	print "  dummy 5"
	pass

def dummy_6args(a,b,c,d,e,f):
	"dummy function, 6 arguments"
	print "  dummy 6"
	pass


class HMM_base:
	"""Generic base class for the discrete ghmm models. """
	def __init__(self):
		"Constructor is to be defined in derived classes."
		pass
		
	def	generate_sequence(self,rnd_seed,length_seq,num_seq, func = dummy_5args):
		"Generate 'num_seq' sequences with length 'length_seq'."
		return sequence(func(self.model,rnd_seed,length_seq,num_seq,100),self.lisalph)

	def make_viterbi(self,mysequence, func = dummy_4args):
		"Takes one sequence and calculates Viterbi-Path and Probability of this Viterbi-Path for this sequence and the HMM."
		log_p = double_array(1)
		viterbi_path = func(self.model,mysequence.seq,mysequence.length,log_p)
		#l = get_int_array_length(viterbi_path)
		#print l
		viterbi_prob = get_arrayd( log_p, 0 )
		
		clean_path = []
		for i in range(mysequence.length * self.model.N):
 			d = get_arrayint(viterbi_path,i)
  			if d >= 0:
	  			clean_path.append(d)
				
		
		return (clean_path,viterbi_prob)
		
	def make_forward(self,mysequence,func = dummy_6args):
		"""Calculates Forward-Probability 'alpha'(t,i) for the intseq object 'mysequence', 
		returns the log probability."""
		log_p = double_array(1)		
		i = self.model.N
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
		
		
		
	def model_prob_dist(self,model, maxT, symmetric, verbose,func = dummy_5args):
		"""Returns the probabilistic distance between 'self' and 'model'"""
		return func(model, model2, maxT, symmetric, verbose)
	
	def get_likelihood(self,sequence,func = dummy_2args):
		"""Returns the log likelihood of 'sequence' under model 'self.model'. """
		return func(self.model, sequence.seq_c) 
		

############################ silent discrete models ###########################################

def read_model_file(fileName, lisalph):
	if (os.path.isfile(fileName)):
		model = graphmlTosdModel(fileName)
		lisalph = lisalph 
		return hmm_model(model, lisalph)
	else :
		raise SystemExit("Error: XML-File "+fileName+" not found!\n") 


class hmm_model(HMM_base):
	"Discrete HMM with silent states"

	def __init__(self,matrans,lisalph,maems=None,pi=None, name="hmm_model" ):
		"'matrans' transition matrix or input XML-file, 'lisalph' list of alphabet, 'maems' emission matrix for each state, 'pi' init prob"
        
		if (type(matrans)!=list):
			self.model = matrans
			self.lisalph = lisalph #.insert(0,'X')
			
		else:
			epsilon = 1e-20
			silent_states = [] 
			silent_flag = 0
			self.lisalph = lisalph
			self.model = model()
			self.model.N = len(matrans)
			self.model.M = len(lisalph)
  
			self.model.prior = -1
			self.model.name = name
			states = arraystate(self.model.N)

			#initialize states
			for i in range(self.model.N):
				state = get_stateptr(states,i)
				state.b = list2arrayd(maems[i])
				state.pi = pi[i]
				if (sum_vector(maems[i]) <= epsilon):
					silent_states.append(1)
					silent_flag = 4
				else:
					silent_states.append(0)

				#set out probabilities
				trans = extract_out(matrans[i])
				state.out_states = trans[0]
				state.out_id = trans[1]
				state.out_a = trans[2]
				#set "in" probabilities
				trans = extract_out(get_col(matrans,i))
				state.in_states = trans[0]
				state.in_id = trans[1]
				state.in_a = trans[2]
				#fix probabilities by reestimation, else 0
				state.fix = 0

			#append states to model
			self.model.s = states
			#self.alloc_silent
			self.model.silent = int_array(self.model.N)
			self.model.model_type = silent_flag
			plist2intarray(self.model.silent, silent_states, self.model.N)
            

	def __str__(self):
		"defines string representation"
		hmm = self.model
		strout = "\nOverview of HMM:\n"
		strout += "\nNumber of states: "+ str(hmm.N)
		strout += "\nSize of Alphabet: "+ str(hmm.M)
		strout += "\nAlphabet: "+str(self.lisalph)

		for k in range(hmm.N):
			state = get_stateptr(hmm.s, k)
			strout += "\n\nState number "+ str(k) +":"
			strout += "\nInitial probability: " + str(state.pi)
			strout += "\nsilent state: " + str(get_arrayint(self.model.silent,k))
			strout += "\nint fix:" + str(state.fix) + "\n"
			strout += "\nOutput probabilites: "
			for outp in range(hmm.M):
				strout+=str(get_arrayd(state.b,outp))+", "
			strout += "\nOutgoing transitions:"
			for i in range( state.out_states):
				strout += "\n  transition to node " + str( get_arrayint(state.out_id,i) ) + " with probability " + str(get_arrayd(state.out_a,i))
			strout +=  "\n\nIngoing transitions:"
			for i in range(state.in_states):
				strout +=  "\n  transition from node " + str( get_arrayint(state.in_id,i) ) + " with probability " + str(get_arrayd(state.in_a,i))
				
		strout += "\n\nSilent states: \n"
		for k in range(hmm.N):
			strout += str(get_arrayint(self.model.silent,k)) + ", "
		strout += "\n"
		return strout
		 
	# ghmm function calls extending the definitons in HMM_base with the apporpriate function names
	def generate_sequence(self,rnd_seed,length_seq,num_seq):
		"generate <num_seq> sequences with length <length_seq>"
		"print generate sequence"
		return HMM_base.generate_sequence(self,rnd_seed,length_seq,num_seq,model_generate_sequences)
	
	def make_viterbi(self,mysequence):
		"takes one sequence and calculates Viterbi-Path and Probability of this Viterbi-Path for this sequence and the HMM"
		return HMM_base.make_viterbi(self, mysequence,viterbi)
	
	def make_forward(self,mysequence):
		"""Calculates Forward-Probability 'alpha'(t,i) for the intseq object 'mysequence', 
		returns the log probability."""
		return HMM_base.make_forward(self, mysequence, foba_forward)
	
	def model_prob_dist(self,model, maxT, symmetric, verbose):
		"""Returns the probabilistic distance between 'self' and 'model'"""
		return HMM_base.model_prob_dist(self,model, maxT, symmetric, verbose, model_prob_dist)
	
	def train_model(self,sequence):
		""" Trains 'self.model' with the sequences contained in sequence object 'sequences'
		 using  reestimate_baum_welch (reestimate.c)."""
		reestimate_baum_welch(self.model, sequence.seq_c)
	
	
	# Equivalent to make_forward
	#def get_likelihood(self,sequence):
	#	"""Returns the log likelihood of 'sequence' under model 'self.model'. """
	#	return HMM_base.get_likelihood(self, sequence, model_likelihood)
		
############################ silent discrete models with transition classes ############################################

def read_sdmodel_file(fileName, lisalph):
	"Read a model file in GXML format."
	if (os.path.isfile(fileName)):
		model = graphmlTosdModel(fileName)
		lisalph = lisalph 
		return hmm_sdmodel(model, lisalph)
	else :
		raise SystemExit("Error: XML-File "+fileName+" not found!\n") 


class hmm_sdmodel(HMM_base):
	"discrete HMM with transitions classes, Wrapper for sdmodel.c"

	def __init__(self,matrans,lisalph,maems=None,pi=None ):
		"<matrans> transition matrix or input XML-file,<lisalph> list of alphabet,<cos> number of output classes,<maems> emission matrix for each state, <pi> init prob"

		if (type(matrans)==str): # if file was specified
			self.model = graphmlTosdModel(matrans)
			self.lisalph = lisalph
			setSwitchingFunction(self.model)
				
		else:
			self.lisalph = lisalph
			self.model = sdmodel()
			self.model.N = len(matrans[0])
			self.model.M = len(lisalph)
			self.model.prior = -1
			cos = len(matrans)
			self.model.cos=cos  # number of classes in GHMM
			# assigning cp_class_change from ghmmwrapper.i as switching function
			setSwitchingFunction( self.model )
			states = arraysdstate(self.model.N)
			#code for silent stuff - necessary
			self.model.silent = int_array(self.model.N)
			bSilent = 0 # exists at least one silent state?

			#initialize states
			for i in range(self.model.N):
				state = get_sdstateptr(states,i)
				state.b = list2arrayd(maems[i])
				state.pi = pi[i]

				#set out probabilities
				tmp = []					
				for j in range(cos):
					tmp.append(matrans[j][i])
				trans = extract_out_probs(tmp,self.model.cos)
				state.out_states = trans[0]
				state.out_id = trans[1]
				state.out_a = trans[2].array

				#set "in" probabilities
				trans = extract_out_probs(get_3d_col(matrans,i),self.model.cos)
				state.in_states = trans[0]
				state.in_id = trans[1]
				state.in_a = trans[2].array
			
				#fix probabilities by reestimation, else 0
				state.fix = 0
				#set if silent
				set_arrayint(self.model.silent,i,self.check_silent(maems[i]))
				bSilent = bSilent or get_arrayint(self.model.silent,i)

			#does model include silent states
			if bSilent:
				self.model.model_type = kSilentStates
				print "Model with silent states"
			else:
				self.model.model_type = kNotSpecified

			#append states to model
			self.model.s = states

	def __str__(self):
		"defines string representation"
		hmm = self.model
		strout = "\nOverview of HMM:"
		strout += "\nNumber of states: "+ str(hmm.N)
		strout += "\nSize of Alphabet: "+ str(hmm.M)
		strout += "\nAlphabet: "+str(self.lisalph)
		strout += "\nNumber of Output-Classes: "+str(hmm.cos)
		print "basic params obtained!"

		for k in range(hmm.N):
			state = get_sdstateptr(hmm.s, k)
			print "accessing state "+str(k)
			strout += "\n\nState number "+ str(k) +":"
			strout += "\nInitial probability: " + str(state.pi)
			strout += "\nOutput probabilites: "
			for outp in range(hmm.M):
				strout+=str(get_arrayd(state.b,outp))+", "
			print "obtained its output probabilites!"
			strout += "\nOutgoing transitions:"
			for i in range( state.out_states):
				strout += "\ntransition to node " + str( get_arrayint(state.out_id,i) )
				for j in range(hmm.cos):
					 strout += "\n\tin class "+str(j)+" with probablity = "+ str(get_2d_arrayd(state.out_a,j,i))
			print "obtained its outgoing transitions!"
			strout +=  "\nIngoing transitions:"
			for i in range(state.in_states):
				strout += "\ntransition from node " + str( get_arrayint(state.in_id,i) )
				for j in range(hmm.cos):
					 print "accessing state.in["+str(i)+","+str(j)+"]"
					 strout += "\n\tin class "+str(j)+" with probablity = "+ str(get_2d_arrayd(state.in_a,j,i))
			print "obtained its ingoing transitions!"
			strout += "\nint fix:" + str(state.fix) + "\n"
			print "obtained its fix!"
		return strout
	
	
	
	def check_silent(self,vecems):
		"Checks whether the emmission vector 'vecems' is silent."
		b = 1;
		for i in range(len(vecems)):
			b = b and (vecems[i]==0.0)
		return b

	def free(self):
		"Wrapps sdmodel_free in sdmodel.c."
		call_sdmodel_free(self.model)

	# ghmm function calls extending the definitons in HMM_base with the apporpriate function names
	def generate_sequence(self,rnd_seed,length_seq,num_seq):
		"generate <num_seq> sequences with length <length_seq>"
		return HMM_base.generate_sequence(self,rnd_seed,length_seq,num_seq,sdmodel_generate_sequences)
	
	def make_viterbi(self,mysequence):
		"takes one sequence and calculates Viterbi-Path and Probability of this Viterbi-Path for this sequence and the HMM"
		return HMM_base.make_viterbi(self, mysequence,sdviterbi)
		
	def make_forward(self,mysequence):
		"calculates Forward-Probability alpha(t,i) for each Emission mysequence[t] and State i, returns C-Double-Matrix"
		return HMM_base.make_forward(self, mysequence, sdfoba_forward)

	def sdmodel_prob_dist(self,model, maxT, symmetric, verbose):
		"""Returns the probabilistic distance between 'self' and 'model'."""	
		return model_prob_distance(model, model2, maxT, symmetric, verbose)
	
	#def get_likelihood(self,sequence):
	#	"""Returns the log likelihood of 'sequence' under model 'self.model'. """
	#	return HMM_base.get_likelihood(model *mo, sequence_t *sq, sdmodel_likelihood)	
		
		
################################ Continous model ##############################################################


def read_smodel_file(fileName,smod_nr):
	""" reads the first <smod_nr> models from <fileName> and returns a python list of smodel pointers"""
	smodel_list = []
	nrModels = int_array(1)
	model_set = smodel_read(fileName, nrModels)
	
	for i in range(smod_nr):
		smodel_list.append( hmm_smodel(get_smodel_ptr(model_set,i)))
	
	# free nr_models
	
	return smodel_list



class hmm_smodel:
	"""Python wrapper for ghmm shmm
		Parameters: 
			matrans: either python matrix of transition probabilies OR
			pointer to completely allocated and initilized smodel
			mean: vector of mean values for the emission pdfs
			var: vector of variances for the emission pdfs
			pi: initial probabilities of the states 	
			weights: emission pdf weights
	
	"""
	def __init__(self,matrans,mean=None,var = None, pi=None, weigths=None ):
		"Constructor"
		if (type(matrans) != list):

			self.model = matrans
		
		else:
			self.model = smodel()
			self.model.N = len(matrans[0])
			self.model.M = len(mean[0])
			self.model.prior = -1
			cos = len(matrans)
			self.model.cos = cos  # number of classes in GHMM
			states = arraysstate(self.model.N)

			# XXX ? switching function ? XXX
			
			#initialize states
			for i in range(self.model.N):
				state = get_sstate_ptr(states,i)
				state.pi = pi[i]
				
				# allocate arrays of emmission parameters
				state.c = list2arrayd(weigths[i])
				state.mue = list2arrayd(mean[i])
				state.u = list2arrayd(var[i])

				#set out probabilities
				tmp = []					
				for j in range(cos):
					tmp.append(matrans[j][i])
				trans = extract_out_probs(tmp,self.model.cos)
				state.out_states = trans[0]
				state.out_id = trans[1]
				state.out_a = trans[2].array

				#set "in" probabilities
				trans = extract_out_probs(get_3d_col(matrans,i),self.model.cos)
				state.in_states = trans[0]
				state.in_id = trans[1]
				state.in_a = trans[2].array
					
				#fix probabilities by reestimation, else 0
				state.fix = 0

			#append states to model
			self.model.s = states
			
			
	def __str__(self):
		"defines string representation"		
		hmm = self.model
		strout = "\nOverview of HMM:"
		strout += "\nNumber of states: "+ str(hmm.N)
		strout += "\nNumber of output distributions per state: "+ str(hmm.M)
		strout += "\nNumber of transition classes: "+str(hmm.cos)
	
		for k in range(hmm.N):
			state = get_sstate(hmm, k)
			strout += "\n\nState number "+ str(k) +":"
			strout += "\nInitial probability: " + str(state.pi)
			strout += "\n"+ str(hmm.M) + " density function(s):\n"
			
			weight = ""
			mue = ""
			u =  ""

			for outp in range(hmm.M):
				weight += str(get_arrayd(state.c,outp))+", "
				mue += str(get_arrayd(state.mue,outp))+", "
				u += str(get_arrayd(state.u,outp))+", "
			
			strout += "  pdf component weights : " + str(weight) + "\n"
			strout += "  mean vector: " + str(mue) + "\n"
			strout += "  variance vector: " + str(u) + "\n"
			strout += "\nOutgoing transitions:"
			
			for i in range( state.out_states):
				strout += "\ntransition to node " + str( get_arrayint(state.out_id,i) )
				for j in range(hmm.cos):
					 strout += "\n\tin class "+str(j)+" with probablity = "+ str(get_2d_arrayd(state.out_a,j,i))
			strout +=  "\nIngoing transitions:"

			for i in range(state.in_states):
				strout += "\ntransition from node " + str( get_arrayint(state.in_id,i) )
				for j in range(hmm.cos):
					 strout += "\n\tin class "+str(j)+" with probablity = "+ str(get_2d_arrayd(state.in_a,j,i))
			strout += "\nint fix:" + str(state.fix) + "\n"
		return strout
	
	
			
	def free(self):
		"Deallocation routine."
		call_smodel_free(self.model)
		
	def generate_sequence(self, seed, seq_len, seq_number, label, Tmax):
		"Generate a sequence_d object from 'self.model'."
		return sequence_d(smodel_generate_sequences(self.model, seed,seq_len,seq_number,label,Tmax))
	
	def make_viterbi(self,mysequence):
		"takes one sequence and calculates Viterbi-Path and Probability of this sequence."
		log_p = double_array(1)
		viterbi_path = sviterbi(self.model,mysequence.seq,mysequence.length,log_p)
		#l = get_int_array_length(viterbi_path)
		#print l
		viterbi_prob = get_arrayd( log_p, 0 )
		#return (intseq(viterbi_path,mysequence.length, 'Viterbi-Path'),viterbi_prob) ### ?

		# XXX richtige laenge des Viterbi Pfades fehlt noch XXX
		return (doubleseq(viterbi_path, 10, 'Viterbi path'),viterbi_prob)
	
	def make_forward(self,mysequence):
		"""calculates Forward-Probability alpha(t,i) for the sequence_d object mysequence, 
		returns the log probability"""


		log_p = double_array(1)		
		i = self.model.N
		logp_sum = 0
		for seq_nr in range(mysequence.seq_c.seq_number):
			t = get_arrayint(mysequence.seq_c.seq_len,seq_nr)
				
			seq = mysequence[seq_nr]
			error = sfoba_logp(self.model, seq.seq, seq.length,log_p)
			if error == -1:
				print "ERROR: Forward finished with -1: Sequence " + str(seq_nr) + " cannot be build."
				exit
			logp_sum  += get_arrayd(log_p,0)
			#print "Seq " + str(seq_nr) + ": " + str(get_arrayd(log_p,0))
			
		return logp_sum # XXX TEST

	
	def Score(self, sequenceSet):
		""" Caller owns this likelihood array """
		print "seqs in sequence set", sequenceSet.seq_number
		self.likelihood = double_array(sequenceSet.seq_number)

		
		result = smodel_individual_likelihoods(self.model, sequenceSet, self.likelihood)
		#if result <= 0:
		#	print "SHMM.Score smodel_individual_likelihoods returns", result
		#	return self.likelihood
		#else:
		return self.likelihood

	def getCState(self, k):
		""" Return a pointer to model->s[k] (sstate)"""
		if ( k >= 0 and k < self.model.N):
			return get_sstate(self.model, k)
		else:
			raise("IndexError")

	def getMean(self, k):
		""" Return a vector of Means"""
		muev = [0.0] * self.model.M
		for k in range(self.model.M):
			muev[k]= get_arrayd(self.model.s[k].mue,k)
		return muev

	def getVariance(self, k):
		""" Return a vector of variance"""
		v = [0.0] * self.model.M
		for k in range(self.model.M):
			sstate = get_sstate(self.model,k)
			v[k]= get_arrayd(sstate.u,k)
		return v
	
	def setMean(self, k, cptMean):
		"""  Mean must be a double array , k is the index """
		smodel_set_mean(self.model, k, cptMean)

	def setVariance(self, k, cptVar):
		""" Variance Mean must be a double array, k is the index  """
		smodel_set_variance(self.model, k, cptVar)

	def setTransition(self, i, j, cos, prob):
		smodel_set_transition(self.model, i,j, cos, float(prob))

	def setPiVector(self, k, pi): 
		smodel_set_pivector(self.model,k,float(pi))

	def setFixVector(self, k, fixv): 
		smodel_set_fixvector(self.model,k,float(fixv))

	def write_hmm_smodel(self, filename):
		""" Write the model to a file in text format """
		call_smodel_print(filename,self.model)
				
# smodel_prob_distance for two smodel pointers


def prob_distance_smodel(model1, model2,seq_length):
	" static wrapper for smodel_prob_distance	"
	return smodel_prob_distance(model1,model2,seq_length,0,0)
		
				
				
################################ Clustering of  smodels ##############################################################

class smodel_cluster_t:
	"Wrapper for the scluster_t struct in scluster.c."	
	def __init__(self,smodFile,max_model=0):
		"Constructor"	
		
		self.scluster_t = scluster_t()
		# reading smodels from file
		mod_nr = int_array(1)
		self.scluster_t.smo = smodel_read(smodFile,mod_nr)
		models_in_file = get_arrayint(mod_nr,0)
		free_arrayi(mod_nr)
		self.scluster_t.smo_number = models_in_file
		
		if models_in_file > max_model:
			print "Warning: models in file > max_models."
			max_nr = models_in_file
		else:
			max_nr = max_model
			smo_array = smodel_array(max_model)
	
			# copying smodel array into larger array to allow modfication during clustering
			for i in range(models_in_file):
				set_smodel_ptr(smo_array,get_smodel_ptr(self.scluster_t.smo,i),i)
	
			#smodel_free(self.sclustering.scluster_t.smo)
			self.scluster_t.smo = smo_array	
						
		self.scluster_t.smo_seq = sequence_d_t_array(max_nr)
		
		self.scluster_t.seq_counter = long_array(max_nr)
		self.scluster_t.smo_Z_MD = double_array(max_nr)
		self.scluster_t.smo_Z_MAW = double_array(max_nr)

		for i in range(max_nr):
			set_arrayd(self.scluster_t.smo_Z_MD,i, 0.0)
			set_arrayd(self.scluster_t.smo_Z_MAW,i, 0.0)
			set_arrayl(self.scluster_t.seq_counter,i,0)

	def __str__(self):
		"defines string representation"
		string = "\nscluster_t:\n"
		string += "Number of models: " + str(self.scluster_t.smo_number) + "\n"
		
		for i in range(self.scluster_t.smo_number):
			seq = get_seq_d_ptr(self.scluster_t.smo_seq,i)
			pseq = sequence_d(seq)
			string += str(pseq) + "\n"
		
		return string
			
			
			
		
class reestimate_smosqd_t:
	""" Wrapper for the smosqd_t struct in sreestimate.c. """
	
	def __init__(self):
		"Constructor"	
		self.smosqd_t = smosqd_t()
	
	# Assesor functions
	def set_model(self, smo):
		self.smosqd_t.smo = smo		
	def set_sqd(self, seq):
		self.smosqd_t.sqd = seq
	def set_logp(self, logp):
		self.smosqd_t.logp = logp	
	def set_eps(self, eps):
		self.smosqd_t.eps = eps		
	def set_max_iter(self, max_iter):
		self.smosqd_t.max_iter = max_iter
	
	def get_model(self):
		return self.smosqd_t.smo	
	def get_sqd(self):
		return self.smosqd_t.sqd
	def get_logp(self):
		return self.smosqd_t.logp 	
	def get_eps(self):
		return self.smosqd_t.eps 		
	def get_max_iter(self):
		return self.smosqd_t.max_iter			


def scluster_print(cl,seq,file_names):	
	"wrapps print_scluster "	
	ch_array = char_array(4)
	for i in range(4):
		set_arraychar(ch_array,i,file_names[i])
	print_scluster(cl,seq,ch_array)


class scluster_env:
	"Wrapper for scluster.c."
	def __init__(self, smodFile,seqFile, rand_labels = 0):
		"Constructor"
		
		# maximum values for models and sequence in the clustering are needed for edit operations
		# like add / remove model etc.
		self.MAX_MODELS = 40
		self.MAX_SEQUENCES = 7000
		
		print "*** Note:"
		print "MAX_MODELS = " + str(self.MAX_MODELS)
		print "MAX_SEQUENCES = " + str(self.MAX_SEQUENCES)
		print "Be sure not to exceed these constants !"
		
		bw_eps = 0.0
		bw_max_iter = 20
		self.sclustering = smodel_cluster_t(smodFile, self.MAX_MODELS)
		
		assert self.sclustering.scluster_t.smo_number <= self.MAX_MODELS
		
		# counter for the maximum number of models present at any time during the clustering,
		# used to assign fresh IDs to new models
		self.most_models = self.sclustering.scluster_t.smo_number
		
		self.seq = seq_d_read(seqFile) # wrapper function seq_d_read calls sequence_d_read
		
		assert self.seq.seq_number <= self.MAX_SEQUENCES
				
		if rand_labels == 1:
			scluster_random_labels(self.seq, self.sclustering.scluster_t.smo_number) ### C-Function call 
			
		
		#print sclustering.scluster_t.smo_number, seq.seq_number
		self.all_log_p = double_2d_array(self.MAX_MODELS,self.MAX_SEQUENCES)
	
		# initializing smosqd_t struct
		self.smosqd_t = reestimate_smosqd_t()
		self.smosqd_t.set_eps(bw_eps)
		self.smosqd_t.set_max_iter(bw_max_iter)
		
		# initializing index map for models
		self.index_map = {} # mapping index in smodel array -> unique smodel ID 
		self.ID_map = {}  #  mapping unique smodel ID -> index in smodel array
		for i in range(self.sclustering.scluster_t.smo_number):
			self.index_map[i] = i
			self.ID_map[i] = i

	def cluster_onestep(self,priors_fixed=0):
		""" One itertation of scluster.c clustering,
		 Arguments:
	 		init_res: an instance of the list returned by smodel_clustering_init. 
			priors_fixed: flag for prior reestimation, 1 = yea, 0 = nay
		 """

		sclustering = self.sclustering
		seq = self.seq
		reest_struct = self.smosqd_t
		all_log_p = self.all_log_p


		changes = 1
		best_logp = double_array(1) 

		# reset termination flag
		changes = 0
		
		print "\n*** start iteration ***\n"
		for i in range(sclustering.scluster_t.smo_number):
			# assigning models to smosqd_t struct
			smodel = get_smodel_ptr(sclustering.scluster_t.smo,i)
			
			print "  model " +str(self.index_map[i] ) + " prior " + str(smodel.prior)
			
			reest_struct.set_model(smodel)
			reest_struct.set_sqd(seq)
			# one row of all_logp is assigned to smosqd_t.logp field
			reest_struct.set_logp(get_row_pointer_d(all_log_p, i))
			
			# calculating the log probabilites for the current model
			scluster_prob(reest_struct.smosqd_t) ### C-Function call
			
			# resetting seq_counter for next iteration
			set_arrayl(sclustering.scluster_t.seq_counter,i,0)

		print 
		double_2d_print(all_log_p,sclustering.scluster_t.smo_number,seq.seq_number)
		print 
		
		# assigning cluster labels to sequences
		for j in range(seq.seq_number):
			label = scluster_best_model(sclustering.scluster_t,j, all_log_p,best_logp ) ### C-Function call
			if changes == 0 and label != get_sequence_d_label(seq,j):
				#print str(label) +" to "+ str(get_sequence_d_label(seq,j))
				changes = 1
				
			print "   sequence " + str(j) +" assigned to model "+ str(self.index_map[label]) + " (l = "+ str(get_arrayd(best_logp,0)) + ")"
			set_sequence_d_label(seq,j,label)
			# cl.seq_counter[sqd->seq_label[j]]++; 
			
			# counting sequences assinged to model
			s_nr = get_arrayl(sclustering.scluster_t.seq_counter,label)
			#print label, s_nr+1
			set_arrayl(sclustering.scluster_t.seq_counter,int(label),int(s_nr +1))

			# Update Z_MD
			old_MD = get_arrayd(sclustering.scluster_t.smo_Z_MD,label)
			new_MD = old_MD + (get_arrayd(seq.seq_w,j) * get_2d_arrayd(all_log_p,label,j)  )
			set_arrayd(sclustering.scluster_t.smo_Z_MD,label, new_MD)
			
			# Update Z_MAW	
			log_apo = double_array(1)
			scluster_log_aposteriori(sclustering.scluster_t,seq,j,log_apo)
			old_MAW = get_arrayd(sclustering.scluster_t.smo_Z_MAW,label)
			new_MAW = old_MAW +  (get_arrayd(seq.seq_w,j) * get_arrayd(log_apo,0 ))
			set_arrayd(sclustering.scluster_t.smo_Z_MAW,label,new_MAW)
		
		# avoid empty models
		scluster_avoid_empty_smodel(seq, sclustering.scluster_t) ### C-Function call
		
		scluster_update(sclustering.scluster_t, seq) ### C-Function call

		for i in range(sclustering.scluster_t.smo_number):
			# assigning models to smosqd_t struct
			smodel = get_smodel_ptr(sclustering.scluster_t.smo,i)
			reest_struct.set_model(smodel)
			aseq = get_seq_d_ptr(sclustering.scluster_t.smo_seq, i)
			reest_struct.set_sqd(seq)
			# one row of all_logp is assigned to smosqd_t.logp field
			reest_struct.set_logp(get_row_pointer_d(all_log_p, i))
		
			
			print "   reestimating model " + str(self.index_map[i])
			# model reestimation
			rs = sreestimate_baum_welch(reest_struct.smosqd_t) ### C-Function call
			if (rs == -1):
				print "Warning: Baum Welch terminated with " + str(rs)
			#print "\nnachher:\n"
			#smodel_print_stdout(reest_struct.smosqd_t.smo)
		
			# updating prior probabilites   			
			if priors_fixed == 0:
				reest_struct.smosqd_t.smo.prior = aseq.total_w / seq.total_w
				#print "prior = " + str(reest_struct.smosqd_t.smo.prior)
			
			# assigning updated model to scluster struct
			set_smodel_ptr(sclustering.scluster_t.smo,reest_struct.smosqd_t.smo,i)
				
		# freeing stuff
		free_arrayd(best_logp)
		
		return changes

	def get_model(self,smo_id):
		"Returns a pointer to the smodel with index 'smo_id'."
		index = self.ID_map[smo_id] 
		return get_smodel_ptr(self.sclustering.scluster_t.smo,index)

	def free_scluster(self):
		"Deallocation."		
		print "freeing sequence_d_t struct seq"
		## Freeing sequence_d_t struct
		#print self.seq
		free_sequence_d(self.seq)
		
		## Freeing scluster_t struct
		# freeing models
		for i in range(self.sclustering.scluster_t.smo_number):
			smodel = get_smodel_ptr(self.sclustering.scluster_t.smo,i)
			call_smodel_free(smodel)
		
		free_smodel_array(self.sclustering.scluster_t.smo)
				
		# freeing sequences
		#sequence_d_free(self.sclustering.scluster_t.smo_seq)		 
		# freeing sequence counters
		free_arrayl(self.sclustering.scluster_t.seq_counter)
		# freeing error function arrays
		free_arrayd(self.sclustering.scluster_t.smo_Z_MD)
		free_arrayd(self.sclustering.scluster_t.smo_Z_MAW)		
		
		## Freeing smosqd_t struct
		free_smosqd_t(self.smosqd_t.smosqd_t)
		
		## Freeing all_logp matrix
		free_2darrayd(self.all_log_p,self.sclustering.scluster_t.smo_number)
		
			
	def add_sequence_d(self, source_seqs):
		"Adds sequences to the clustering.'source_seqs' is a sequence_d object."
		sequence_d_add(self.seq,source_seqs.seq_c)
		
	def add_SHMM_model(self, smodel,prior):
		"Adds a smodel to the clustering. 'smodel' is a pointer to a smodel struct."
		smodel.model.prior = prior
		set_smodel_ptr(self.sclustering.scluster_t.smo, smodel.model, self.sclustering.scluster_t.smo_number ) 
		
		
		# rescaling prior probabilities
		scale = 1/(1-prior)
		for i in range(self.sclustering.scluster_t.smo_number):
			smodel = get_smodel_ptr(self.sclustering.scluster_t.smo,i)
			smodel.prior = smodel.prior / scale
		
		
		print " new model added at index " + str(self.sclustering.scluster_t.smo_number) +" with ID "+ str(self.most_models)
		
		self.index_map[self.sclustering.scluster_t.smo_number] = self.most_models
		self.ID_map[self.most_models] = self.sclustering.scluster_t.smo_number
		
		self.most_models += 1
		
		# incrementing number of models
		new_nr = self.sclustering.scluster_t.smo_number +1
		self.sclustering.scluster_t.smo_number = new_nr

	def remove_model(self, smo_id, rescale_priors = 1):
		"Deletes the smodel with model ID 'smo_id' from the clustering."
		index = self.ID_map[smo_id]
		assert index < self.sclustering.scluster_t.smo_number
		scale = 1 - get_smodel_ptr(self.sclustering.scluster_t.smo,index).prior

		# deleting model		
		call_smodel_free(get_smodel_ptr(self.sclustering.scluster_t.smo,index))
		# assigning former last model to the index we are goign to delete
		last_model = self.sclustering.scluster_t.smo_number -1
		set_smodel_ptr(self.sclustering.scluster_t.smo, get_smodel_ptr(self.sclustering.scluster_t.smo,last_model),index)
		
		# decrementing number of models
		#n = self.sclustering.scluster_t.smo_number
		self.sclustering.scluster_t.smo_number = last_model 
		
		# assigning new model index to index_map
		self.index_map[index] = self.index_map[last_model]
		self.ID_map[self.index_map[last_model]] = index
	
		#print self.index_map[index], self.ID_map[last_model],index
		
		del self.index_map[last_model]
		del self.ID_map[smo_id]
		
		
		
		if rescale_priors == 1:
			# rescaling prior probabilities		
			for i in range(last_model):
				smodel = get_smodel_ptr(self.sclustering.scluster_t.smo,i)
				smodel.prior = smodel.prior / scale

	def merge_models(self, smo_id1, smo_id2):
		"""Merges models with IDs 'smo_id1' and 'smo_id2' by reestimating the model with more sequences
			attached with the sequences of both models.
		"""
		index1 = self.ID_map[smo_id1]
		index2 = self.ID_map[smo_id2]
		
		# get pointers to models to be merged
		smo_num1 = get_arrayl(self.sclustering.scluster_t.seq_counter,index1)
		smo_num2 = get_arrayl(self.sclustering.scluster_t.seq_counter,index2)
	
		# the merging is implemented as the retraining of the model with the greater number of assinged sequences 
		# with the sequences of both models 
		if smo_num1 >= smo_num2:
			model2retrain = get_smodel_ptr(self.sclustering.scluster_t.smo,index1)
			model2dump = get_smodel_ptr(self.sclustering.scluster_t.smo,index2)	
			retrain_index = index1
			dump_index = index2
			dump_id = smo_id2
			retrain_id = smo_id1
		else:
			model2retrain = get_smodel_ptr(self.sclustering.scluster_t.smo,index2)
			model2dump = get_smodel_ptr(self.sclustering.scluster_t.smo,index1)	
			retrain_index = index2
			dump_index = index1
			retrain_id = smo_id2
			dump_id = smo_id1
			
		# get sequence_d_t pointers 
		retrain_seq = get_seq_d_ptr(self.sclustering.scluster_t.smo_seq,retrain_index)	
		dump_seq = get_seq_d_ptr(self.sclustering.scluster_t.smo_seq,dump_index)	
		# adding sequences to model to be retrained
		
		sequence_d_add(retrain_seq,dump_seq)			
		
		# raising sequence counter for this model
		set_arrayl(self.sclustering.scluster_t.seq_counter, retrain_index, int(smo_num1 + smo_num2))

		# inserting model to be retrained into smosqd_t struct
		self.smosqd_t.set_model(model2retrain)
		self.smosqd_t.set_sqd(retrain_seq)
		self.smosqd_t.set_logp(get_row_pointer_d(self.all_log_p, retrain_index))
	
		# reestimating model
		sreestimate_baum_welch(self.smosqd_t.smosqd_t)	
	

	
		# the merged model gets a fresh ID
		print "setting ID of merged model to " + str(self.most_models)
		self.index_map[retrain_index] = self.most_models
		del self.ID_map[retrain_id]
		self.ID_map[self.most_models] = retrain_index
		
		# adding model priors
		model2retrain.prior += model2dump.prior
		
		self.most_models += 1
		self.remove_model(dump_id,0)
		return self.most_models -1

		
		
	def deviate_model(self, smo_id,trans_dev,mean_dev,var_dev,weight_dev,prior = 0):
		"Introduces random devation into the model with model ID 'smo_id'."
		index = self.ID_map[smo_id]
		base_hmm = get_smodel_ptr(self.sclustering.scluster_t.smo, index)
		hmm = smodel_copy(base_hmm)		
		
				
		for k in range(hmm.N):
			state = get_sstate(hmm, k)
			
			#print "*** STATE " + str(k)
			
			w_sum = 0
			for outp in range(hmm.M):
				mean = get_arrayd(state.mue,outp)
				var = get_arrayd(state.u,outp)
				w = get_arrayd(state.c,outp)
				
				#print mean, var, w
				
				mean_rand = mean + uniform(-mean_dev,mean_dev)
				var_rand = abs(var + uniform(-var_dev,var_dev))
				weight_rand = w + uniform(-weight_dev,weight_dev)
				if weight_rand < 0:
					weight_rand = 1e-6
				w_sum += weight_rand
				
				#print mean_rand, var_rand, weight_rand
				
				set_arrayd(state.mue,outp,mean_rand)
				set_arrayd(state.u, outp,var_rand)
				set_arrayd(state.c, outp,weight_rand)
				
			# renormalizing weight vector	
			for outp in range(hmm.M):
				w = get_arrayd(state.c,outp)
				set_arrayd(state.c, outp,w/w_sum)	
			
			#print get_arrayd(state.c,0)
			#print " transitions"
			
			for tclass in range(hmm.cos):
				p_sum = 0
				for sid in range( state.out_states):
					out_id = get_arrayint(state.out_id,sid) 
				
					p = get_2d_arrayd(state.out_a,tclass,sid)

					if p > 0:
						trans_rand = p + uniform(-trans_dev,trans_dev)
						if trans_rand < 0:
							trans_rand = 0
							print "transitions smaller 0 set to 0."
						#print "  ",p," -->  ", trans_rand
						
						set_2d_arrayd(state.out_a,tclass,sid,trans_rand)
						p_sum += trans_rand
			
			
				# renormalizing probabilities
				for sid in range( state.out_states):
					p = get_2d_arrayd(state.out_a,tclass,sid)
					if p_sum != 0:
						set_2d_arrayd(state.out_a,tclass, sid, p/p_sum)
				
				
				# updating the in_a probabilites to the deviated values
				for sid in range( state.out_states):
					out_id = get_arrayint(state.out_id,sid) 
					in_state = get_sstate(hmm, out_id)
					p_new = get_2d_arrayd(state.out_a,tclass,sid)
										
					for t in range(in_state.in_states):
						i = get_arrayint(in_state.in_id,t)
						if get_arrayint(in_state.in_id,t) == k:
							set_2d_arrayd(in_state.in_a,tclass,t,p_new)
							break
		
		m = hmm_smodel("bla",0,hmm,1) 
		self.add_SHMM_model( m,prior)
							
	def update_labels(self):
		"Assigns the model IDs to the corresponding sequences."
		for i in range(self.seq.seq_number):
			old_label = get_arrayl(self.seq.seq_label,i)
			new_label = self.index_map[old_label]
			set_arrayl(self.seq.seq_label,i,new_label)
			
			
			
			
			
			
