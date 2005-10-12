/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/gauss_tail.h
*       Authors:  Bernhard Knab
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


#ifndef GHMM_GAUSS_TAIL_H
#define GHMM_GAUSS_TAIL_H


#ifdef __cplusplus
extern "C" {
#endif

/**
   @name some calculations concerning gaussian tail function
 */

/*@{ */
/**
 */
  double ighmm_gtail_pmue (double mue, double A, double B, double eps);

/**
 */
  double ighmm_gtail_pmue_umin (double mue, double A, double B, double eps);

/** 
    @name Function to find the roots of the truncated normal density function.
*/
  double ighmm_gtail_pmue_interpol (double mue, double A, double B, double eps);

#ifdef __cplusplus
}
#endif
/*@} */
#endif                          /* GHMM_GAUSS_TAIL_H */
