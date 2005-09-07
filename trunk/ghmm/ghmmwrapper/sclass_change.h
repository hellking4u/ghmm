/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: sclass_change.h
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
#include <ghmm/sdmodel.h>

/* Function for class changes for switching models

   smo is a smodel struct
   seq is a double matrix of observations 
   k is the first current sequence index (corresponding to the rows in seq)
   t is the second current sequence index (corresponding to seq[k] )
*/
int cp_class_change( smodel *smo, int *seq, int k, int t);

/*
   setSwitchingFunction assigns cp_class_change as switching function in model smo.
   Needs to be modified for user defined C switching function.
*/
void setSwitchingFunction( smodel *smd );

/* Assignment of Python module and function for class change. The values are stored in smo->class_change.
   
   smo: smodel struct with multiple transition classes
   python_module: Name of the module the switching function is defined in
   python_function: Name of the Python function to be used. 
   IMPORTANT: python_function must have the same signature as cp_class_change (that means three arguments:
   first the sequence, second the sequence index, third the time step in the current sequence. See class_change.py in the
   ghmmwrapper directory for an example.)

*/
void setPythonSwitching( smodel *smd, char* python_module, char* python_function);

/* Implements the Python Callback to the switching function defined in smo->class_change.
   Arguments are identical to cp_class_change (s.a.)

*/
int python_class_change( smodel* smo, int* seq, int k, int t );










