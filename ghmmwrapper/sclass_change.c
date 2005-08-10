/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: sclass_change.c
*       Authors:  Benjamin Georgi
*
*       Copyright (C) 1998-2004 Alexander Schliep
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
*                               Berlin
*
*       Contact: schliep@ghmm.org
*
*       This library is free software; you can redistribute it and/or
*       modify it under the terms of the GNU Library General Public
*       License as published by the Free Software Foundation; either
*       version 2 of the License, or (at your option) any later version.
*
*       This library is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*       Library General Public License for more details.
*
*       You should have received a copy of the GNU Library General Public
*       License along with this library; if not, write to the Free
*       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*
*       This file is version $Revision$
*                       from $Date$
*             last change by $Author$.
*
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/smodel.h>
#include <ghmm/mes.h>
#include <Python.h>

int cp_class_change( smodel *smo, double* seq, int k, int t) {
  int res;
  
  if((t%2) == 0){
    /*if(t < 500){*/
      res= 0;
  }
  else {
      res = 1;
  }        
  printf("cp_class_change with value %d -> return %d\n",t,res);    

  return res;
} 		
		

void setSwitchingFunction( smodel *smd ) {
  smd->class_change->get_class = cp_class_change;
}


int python_class_change( smodel* smo, double *seq, int k, int t ){
   char* ModuleName = smo->class_change->python_module;
   char* FunctionName = smo->class_change->python_function;
   int class,i;

   PyObject *pName,  *pDict, *pArgs, *pValue, *pList;
  
   /* Python module and function are static */
   static PyObject *pModule = NULL;
   static PyObject *pFunc = NULL;
   
   /* importing module and function on first function call */
   if (pModule == NULL) {   
     printf("C: Importing Python module ... ");
     pName = PyString_FromString(ModuleName);
   
     pModule = PyImport_Import(pName);       // Import module
     if(!pModule) {
       printf("python_class_change ERROR: Module %s not found.\n",ModuleName);
       return(-1);
     }    
   
     pDict = PyModule_GetDict(pModule);
     printf("done.\n");    
    
     //printf("C: Calling Python with value %d\n",t);
     pFunc = PyDict_GetItemString(pDict, FunctionName);
     if(!pFunc) {
       printf("python_class_change ERROR: Function %s not found.\n",FunctionName);
       return(-1);
     }  
     Py_DECREF(pDict); 
     Py_DECREF(pName); 
   }
   
   pArgs = PyTuple_New(3);
   
   pList = PyList_New(t);
   for(i=0;i<t;i++){
     pValue = PyFloat_FromDouble(seq[i]);
     PyList_SetItem(pList, i, pValue);
   }    
   PyTuple_SetItem(pArgs, 0, pList); 

   pValue = PyInt_FromLong((long)k);
   PyTuple_SetItem(pArgs, 1, pValue); 

   pValue = PyInt_FromLong((long)t);
   PyTuple_SetItem(pArgs, 2, pValue); 
   
   pValue = PyObject_CallObject(pFunc, pArgs); // Calling Python 

   /* parsing the result from Python to C data type */
   class = PyInt_AsLong(pValue);
  /*printf("C: The returned class is %d\n",class);*/
     
   /* cleaning up */
   Py_DECREF(pArgs); 
   Py_DECREF(pValue); 
   Py_DECREF(pList);    
   
   return class; 
 
}


void setPythonSwitching( smodel *smd, char* python_module, char* python_function ){
   if(!smd->class_change) {
     printf("setPythonSwitching ERROR: class_change struct not initialized.\n");
   }     
   smd->class_change->python_module = python_module;
   smd->class_change->python_function = python_function;     
   smd->class_change->get_class = python_class_change ;
}
