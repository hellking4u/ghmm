/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/discrime.h
*       Authors:  Janne Grunau
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


#ifndef DISCRIME_H
#define DISCRIME_H

#ifdef __cplusplus
extern "C" {
#endif

model** discrime_modelarray_alloc(int size);
void discrime_modelarray_dealloc(model** mos);

void discrime_modelarray_setptr(model** mos, model* mo, int pos);
model* discrime_modelarray_getptr(model** mos, int pos);

sequence_t** discrime_seqarray_alloc(int size);
void discrime_seqarray_dealloc(sequence_t** seqs);

void discrime_seqarray_setptr(sequence_t** seqs, sequence_t* seq, int pos);
sequence_t* discrime_seqarray_getptr(sequence_t** seqs, int pos);

/*----------------------------------------------------------------------------*/
/**
   Trains two or more models to opimise the discrimination between the
   classes in the trainingset.
   The models must have the same topology. (checked)
   @return                 0/-1 success/error
   @param mo:              array of pointers to some models
   @param sqs:             array of annotated sequence sets
   @param noC:             number of classes
   @param gradient:        if gradient == 0 try a closed form solution
                           otherwise a gradient descent
 */
extern int discriminative(model** mo, sequence_t** sqs, int noC, int gradient);

/*----------------------------------------------------------------------------*/
/**
   Returns the value of the in this discriminative training algorithm optimised
   function for a tupel of HMMs and sequencesets.
   @return                 value of funcion
   @param mo:              array of pointers to some models
   @param sqs:             array of annotated sequence sets
   @param noC:             number of classes
*/
double discrime_compute_performance(model** mo, sequence_t** sqs, int noC);

void discrime_print_statistics(model** mo, sequence_t** sqs, int noC,
			       int* falseP, int* falseN);

#ifdef __cplusplus
}
#endif

#endif
