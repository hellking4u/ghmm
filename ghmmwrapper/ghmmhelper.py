import ghmmwrapper
import math
import os.path
from modhmmer import *
from random import *

def emslist2array(lisems):
    "converts emission list to C double array"
    arrems = ghmmwrapper.double_array(len(lisems))
    for i in range(len(lisems)):
        ghmmwrapper.set_arrayd(arrems,i,lisems[i])
    return arrems


def extract_out(lisprobs):
    """ Helper function for building HMMs from matrices: Used for
	    transition matrices without  transition classes.

		Extract out-/ingoing transitions from the row-vector resp. the
		column vector (corresponding to incoming transitions) of the
		transition matrix

		Allocates: .[out|in]_id and .[out|in]_a vectors
    """
    lis = []
    for i in range(len(lisprobs)):
        if lisprobs[i]!=0:
            lis.append(i)
    trans_id = ghmmwrapper.int_array(len(lis))
    trans_prob = ghmmwrapper.double_array(len(lis))
    for i in range(len(lis)):
        ghmmwrapper.set_arrayint(trans_id,i,lis[i])
        ghmmwrapper.set_arrayd(trans_prob,i,lisprobs[lis[i]])
    return [len(lis),trans_id,trans_prob]

def list2arrayd(lisems):
	"converts python list to C double array"
	arrems = ghmmwrapper.double_array(len(lisems))
	for i in range(len(lisems)):
		ghmmwrapper.set_arrayd(arrems,i,lisems[i])
	return arrems


def extract_out_probs(lisprobs,cos):
##    """ Helper function for building HMMs from matrices: Used for
##	    transition matrices with 'cos' transition classes.

##		Extract out-/ingoing transitions from a matric consiting of
##		the row-vectors resp. the column vectors (corresponding to
##		incoming transitions) of the 'cos' transition matrices.
##		Hence, input is a 'cos' x N matrix.

##		Allocates: .[out|in]_id vector and .[out|in]_a array (of size cos x N)
##    """
	lis = []
	# parsing indixes belonging to postive probabilites
	for j in range(cos):
		for i in range(len(lisprobs[0])):
			if lisprobs[j][i] != 0 and i not in lis:
				lis.append(i)
	# print "lis: ", lis			


	trans_id   = ghmmwrapper.int_array(len(lis))
	probsarray = ghmmwrapper.double_2d_array(cos, len(lis)) # C-function
	
	# creating list with positive probabilities
	for k in range(cos):
		for j in range(len(lis)):
			ghmmwrapper.set_2d_arrayd(probsarray,k,j, lisprobs[k][lis[j]])

	trans_prob = twodim_double_array(probsarray, cos, len(lis)) # python CLASS, C internal
	
	#print trans_prob
	
	# initializing c state index array
	for i in range(len(lis)):
		ghmmwrapper.set_arrayint(trans_id,i,lis[i])
	return [len(lis),trans_id,trans_prob]

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
        return ghmmwrapper.get_2d_arrayd(self.array,index[0],index[1])

    def __setitem__(self,index,value):
        "defines twodim_double_array[index[0],index[1]]"
        if (len(index) == 2):
            ghmmwrapper.set_2d_arrayd(self.array,index[0],index[1],value)
                
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


class double_array:
    """A C-double array"""

    def __init__(self, array, columns, columnlabels=None):
        """Constructor"""
        self.array = array
        self.rows = 1
        self.columns = columns
        self.size = columns
        self.columnlabels = columnlabels

    def __getitem__(self,index):
        """defines double_array[index] """
        return ghmmwrapper.get_arrayd(self.array,index)

    def __setitem__(self,index,value):
        """ double_array[index] = value """
        ghmmwrapper.set_arrayd(self.array,index,value)

    def __str__(self):
        """defines string representation"""
        strout = "\n"
        if (self.columnlabels is not None):
            for k in range(len(self.columnlabels)):
                strout+="\t"
                strout+= str(self.columnlabels[k])
                strout += "\n"
        for i in range(self.columns):
            strout += "%2.4f" % self[i]
            strout += "\t"
            strout += "\n"
        return strout

