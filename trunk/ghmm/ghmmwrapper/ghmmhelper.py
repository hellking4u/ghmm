#/*******************************************************************************
#*
#*       This file is part of the General Hidden Markov Model Library,
#*       GHMM version __VERSION__, see http://ghmm.org
#*
#*       Filename: ghmmhelper.py
#*       Authors:  Benjamin Georgi
#*
#*       Copyright (C) 1998-2004 Alexander Schliep
#*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
#*       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
#*                               Berlin
#*
#*       Contact: schliep@ghmm.org
#*
#*       This library is free software; you can redistribute it and/or
#*       modify it under the terms of the GNU Library General Public
#*       License as published by the Free Software Foundation; either
#*       version 2 of the License, or (at your option) any later version.
#*
#*       This library is distributed in the hope that it will be useful,
#*       but WITHOUT ANY WARRANTY; without even the implied warranty of
#*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#*       Library General Public License for more details.
#*
#*       You should have received a copy of the GNU Library General Public
#*       License along with this library; if not, write to the Free
#*       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#*
#*
#*       This file is version $Revision$
#*                       from $Date$
#*             last change by $Author$.
#*
#*******************************************************************************/
import ghmmwrapper
import math
import os.path
from modhmmer import *
from random import *


def list2arrayd(lisems):
    "converts python list to C double array"
    arrems = ghmmwrapper.double_array(len(lisems))
    for i in range(len(lisems)):
        ghmmwrapper.set_arrayd(arrems,i,lisems[i])
    return arrems

def list2arrayint(lisems):
    "converts python list to C int array"
    arrems = ghmmwrapper.int_array(len(lisems))
    for i in range(len(lisems)):
        ghmmwrapper.set_arrayint(arrems,i,lisems[i])
    return arrems


def arrayd2list(carray,length):
    "converts C double array to python list."
    l = []
    for i in range(length):
        l.append(ghmmwrapper.get_arrayd(carray,i))
    return l
    
def arrayint2list(carray,length):
    "converts C double array to python list."
    l = []
    for i in range(length):
        l.append(ghmmwrapper.get_arrayint(carray,i))
    return l
    
def matrixd2list(cmatrix,row,col):
    llist = []
    for i in range(row):
        llist.append([])    
     
    for i in range(row):
        for j in range(col):
            llist[i].append(ghmmwrapper.get_2d_arrayd(cmatrix,i,j) )
    return llist       
        


def list2matrixd(matrix):
    """ Allocation and initialization of a double** based on a   
        two dimensional Python list (list of lists). The number of elements
        in each column can vary.
    """
    rows = len(matrix)
    
    seq = ghmmwrapper.double_2d_array_nocols(rows)
    col_len = []
    for i in range(rows):
        l = len(matrix[i])
        col = ghmmwrapper.double_array(l)
        for j in range(l):
            ghmmwrapper.set_arrayd(col,j,matrix[i][j]) 
        
        ghmmwrapper.set_2d_arrayd_col(seq,i,col)
        col_len.append(l)

    return (seq,col_len)    

def list2matrixint(matrix):
    """ Allocation and initialization of an int** based on a   
        two dimensional Python list (list of lists). The number of elements
        in each column can vary.
    """
    rows = len(matrix)
    
    seq = ghmmwrapper.int_2d_array_nocols(rows)
    col_len = []
    for i in range(rows):
        l = len(matrix[i])
        col = ghmmwrapper.int_array(l)
        for j in range(l):
            ghmmwrapper.set_arrayint(col,j,matrix[i][j]) 
        
        ghmmwrapper.set_2d_arrayint_col(seq,i,col)
        col_len.append(l)

    return (seq,col_len)    

def matrixint2list(cmatrix,row,col):
    llist = []
    for i in range(row):
        llist.append([])    
     
    for i in range(row):
        for j in range(col):
            llist[i].append(ghmmwrapper.get_2d_arrayint(cmatrix,i,j) )
    return llist       
        


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



#def extract_out_probs(lisprobs,cos):
##    """ Helper function for building HMMs from matrices: Used for
##        transition matrices with 'cos' transition classes.
##        Extract out-/ingoing transitions from a matric consiting of
##        the row-vectors resp. the column vectors (corresponding to
##        incoming transitions) of the 'cos' transition matrices.
##        Hence, input is a 'cos' x N matrix.
##        Allocates: .[out|in]_id vector and .[out|in]_a array (of size cos x N)
##    """
#    lis = []
    # parsing indixes belonging to postive probabilites
#    for j in range(cos):
#        for i in range(len(lisprobs[0])):
#            if lisprobs[j][i] != 0 and i not in lis:
#                lis.append(i)
#    print "lis: ", lis            
#    trans_id   = ghmmwrapper.int_array(len(lis))
#    probsarray = ghmmwrapper.double_2d_array(cos, len(lis)) # C-function
    # creating list with positive probabilities
#    for k in range(cos):
#        for j in range(len(lis)):
#            ghmmwrapper.set_2d_arrayd(probsarray,k,j, lisprobs[k][lis[j]])
#    trans_prob = twodim_double_array(probsarray, cos, len(lis)) # python CLASS, C internal
#    
#    print trans_prob
    # initializing c state index array
#    for i in range(len(lis)):
#        ghmmwrapper.set_arrayint(trans_id,i,lis[i])
#    return [len(lis),trans_id,trans_prob]

def extract_out_cos(transmat, cos, state):
    """ Helper function for building HMMs from matrices: Used for
        transition matrices with 'cos' transition classes.

        Extract outgoing transitions for 'state' from the complete list
        of transition matrices

        Allocates: .out_id vector and .out_a array (of size cos x N)
    """
    
    lis = []
    # parsing indixes belonging to postive probabilites
    for j in range(cos):
        for i in range(len(transmat[j][state])):
            if transmat[j][state][i] != 0.0 and i not in lis:
                lis.append(i)

    #lis.sort()
    #print "lis: ", lis            
    
    trans_id   = ghmmwrapper.int_array(len(lis))
    probsarray = ghmmwrapper.double_2d_array(cos, len(lis)) # C-function
    
    # creating list with positive probabilities
    for k in range(cos):
        for j in range(len(lis)):
            ghmmwrapper.set_2d_arrayd(probsarray,k,j, transmat[k][state][lis[j]])
    
    # initializing C state index array
    for i in range(len(lis)):
        ghmmwrapper.set_arrayint(trans_id,i,lis[i])
    return [len(lis),trans_id,probsarray]

def extract_in_cos(transmat, cos, state):
    """ Helper function for building HMMs from matrices: Used for
        transition matrices with 'cos' transition classes.

        Extract ingoing transitions for 'state' from the complete list
        of transition matrices

        Allocates: .in_id vector and .in_a array (of size cos x N)
    """
    lis = []

    # parsing indixes belonging to postive probabilites
    for j in range(cos):
        transmat_col_state = map( lambda x: x[state], transmat[j])
        for i in range(len(transmat_col_state)):
            
            if transmat_col_state[i] != 0.0 and i not in lis:
                lis.append(i)

    #lis.sort()
    #print "lis: ", lis            


    
    trans_id   = ghmmwrapper.int_array(len(lis))
    probsarray = ghmmwrapper.double_2d_array(cos, len(lis)) # C-function
    
    # creating list with positive probabilities
    for k in range(cos):
        for j in range(len(lis)):
            ghmmwrapper.set_2d_arrayd(probsarray,k,j, transmat[k][lis[j]][state])
    
    # initializing C state index array
    for i in range(len(lis)):
        ghmmwrapper.set_arrayint(trans_id,i,lis[i])
    return [len(lis),trans_id,probsarray]

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
        
            
def classNumber(A):
    """ Returns the number of transition classes in the matrix A   """
    cos = 0
    if type(A[0][0]) == list:
        cos = len(A)
    else: 
        cos = 1
    return cos

