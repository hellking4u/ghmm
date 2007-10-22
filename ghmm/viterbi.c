/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/viterbi.c
*       Authors:  Wasinee Rungsarityotin, Benjamin Georgi
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

#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <float.h>
#include <math.h>
#include <assert.h>

#include "ghmm.h"
#include "mes.h"
#include "mprintf.h"
#include "matrix.h"
#include "sdmodel.h"
#include "ghmm_internals.h"

typedef struct local_store_t {
    double **log_in_a;
    double **log_b;
    double *phi;
    double *phi_new;
    int **psi;

    int *topo_order;
    int topo_order_length;
} local_store_t;


/*----------------------------------------------------------------------------*/
static int viterbi_free(local_store_t **v, int n, int len)
{
#define CUR_PROC "viterbi_free"
    int j;
    mes_check_ptr(v, return (-1));

    if (!*v)
        return (0);

    for (j = 0; j < n; j++)
        m_free((*v)->log_in_a[j]);

    m_free((*v)->log_in_a);
    ighmm_cmatrix_free(&((*v)->log_b), n);
    m_free((*v)->phi);
    m_free((*v)->phi_new);
    /*ighmm_dmatrix_free( &((*v)->psi), len ); */
    ighmm_dmatrix_stat_free(&((*v)->psi));
    m_free((*v)->topo_order);
    m_free(*v);

    return (0);
#undef CUR_PROC
}                               /* viterbi_free */

/*----------------------------------------------------------------------------*/
static local_store_t *viterbi_alloc(ghmm_dmodel *mo, int len)
{
#define CUR_PROC "sdviterbi_alloc"
    local_store_t *v = NULL;
    int j;
    ARRAY_CALLOC(v, 1);

    /* Allocate the log_in_a's -> individal lenghts */
    ARRAY_CALLOC(v->log_in_a, mo->N);
    for (j = 0; j < mo->N; j++) {
        ARRAY_CALLOC(v->log_in_a[j], mo->s[j].in_states);
    }

    v->log_b = ighmm_cmatrix_alloc(mo->N, len);
    if (!(v->log_b)) {
        GHMM_LOG_QUEUED(LCONVERTED);
        goto STOP;
    }
    ARRAY_CALLOC(v->phi, mo->N);
    ARRAY_CALLOC(v->phi_new, mo->N);
    v->psi = ighmm_dmatrix_stat_alloc(len, mo->N);
    if (!(v->psi)) {
        GHMM_LOG_QUEUED(LCONVERTED);
        goto STOP;
    }

    v->topo_order_length = 0;
    ARRAY_CALLOC(v->topo_order, mo->N);

    return (v);
  STOP:                        /* Label STOP from ARRAY_[CM]ALLOC */
    viterbi_free(&v, mo->N, len);
    return (NULL);
#undef CUR_PROC
}                               /* viterbi_alloc */

/*----------------------------------------------------------------------------*/
static void Viterbi_precompute(ghmm_dmodel *mo, int *o, int len, local_store_t *v)
{
#define CUR_PROC "viterbi_precompute"
    int i, j, t;

    /* Precomputing the log(a_ij) */
    for (j = 0; j < mo->N; j++)
        for (i = 0; i < mo->s[j].in_states; i++)
            if (mo->s[j].in_a[i] == 0.0)        /* DBL_EPSILON ? */
                v->log_in_a[j][i] = +1; /* Not used any further in the calculations */
            else
                v->log_in_a[j][i] = log(mo->s[j].in_a[i]);


    /* Precomputing the log(bj(ot)) */
    for (j = 0; j < mo->N; j++)
        for (t = 0; t < len; t++) {
            if (o[t] != mo->M) {
                if (mo->s[j].b[o[t]] == 0.0)    /* DBL_EPSILON ? */
                    v->log_b[j][t] = +1;
                else
                    v->log_b[j][t] = log(mo->s[j].b[o[t]]);
            } else {
                v->log_b[j][t] = 0.0;
            }
        }

#undef CUR_PROC
}                               /* viterbi_precompute */

/*----------------------------------------------------------------------------*/
static void __viterbi_silent(ghmm_dmodel *mo, int t, local_store_t *v)
{
#define CUR_PROC "viterbi_silent"
    int topocount;
    int i, k;
    double max_value, value;

    for (topocount = 0; topocount < mo->topo_order_length; topocount++) {
        k = mo->topo_order[topocount];
        if (mo->silent[k]) {    /* Silent states */
            /* Determine the maximum */
            /* max_phi = phi[i] + log_in_a[j][i] ... */
            max_value = -DBL_MAX;
            v->psi[t][k] = -1;
            for (i = 0; i < mo->s[k].in_states; i++) {

                if (v->phi[mo->s[k].in_id[i]] != +1 && v->log_in_a[k][i] != +1) {
                    value = v->phi[mo->s[k].in_id[i]] + v->log_in_a[k][i];
                    if (value > max_value) {
                        max_value = value;
                        v->psi[t][k] = mo->s[k].in_id[i];
                    }
                }
            }
            /* No maximum found (that is, state never reached)
               or the output O[t] = 0.0: */
            if (max_value == -DBL_MAX) {
                v->phi[k] = +1;
            } else {
                v->phi[k] = max_value;
            }
        }
    }
#undef CUR_PROC
}

/*============================================================================*/
/* auxilary function. Reallocates an integer array to a given length
   and initialises the new positions with -1 */
static int extend_int_array(int *array, int cur_len, int extend)
{
#define CUR_PROC "extend_int_array"
    int j;
    cur_len += extend;
    ARRAY_REALLOC(array, cur_len);
    /* initialising new memory with -1  */
    for (j = cur_len - 1; j >= (cur_len - extend); j--) {
        array[j] = -1;
    }
    return (cur_len);
  STOP:
    return (-1);
#undef CUR_PROC
}

/*============================================================================*/
/** Return the viterbi path of the sequence.  */
int *ghmm_dmodel_viterbi(ghmm_dmodel * mo, int *o, int len, double *log_p)
{
#define CUR_PROC "ghmm_dmodel_viterbi"

    int *state_seq = NULL;
    int t, j, i, k, St;
    double value, max_value;
    local_store_t *v;

    /* length_factor determines the size of the memory allocation for the state
       path array in case the model contains silent states.
       The larger length_factor is chosen, the less reallocs will be necessary
       but it increases the amount of allocated memory that is not used. */
    int length_factor = 2;

    /* the lenght of the viterbi path is unknown for models with silent states.
       The maximum length is given by mo->N * len. In order to keep
       memory consumption in check we first allocate lenght_factor * len
       and realloc more space as needed. */
    int len_path;

    int cur_len_path = 0;       /* the current length of the viterbi path */

    int lastemState;
    int state_seq_index;

    /* for silent states: initializing path length with a multiple
       of the sequence length
       and sort the silent states topological */
    if (mo->model_type & GHMM_kSilentStates) {
        len_path = length_factor * len;
        ghmm_dmodel_order_topological(mo);
    }
    /* if there are no silent states, path and sequence length are identical */
    else
        len_path = len;

    /* Allocate the matrices log_in_a, log_b,Vektor phi, phi_new, Matrix psi */
    v = viterbi_alloc(mo, len);
    if (!v) {
        GHMM_LOG_QUEUED(LCONVERTED);
        goto STOP;
    }

    /* allocating state_seq array */
    ARRAY_CALLOC(state_seq, len_path);

    /* initialization of state_seq with -1, only necessary for silent state models */
    /* XXX use memset?, not faster, assumes 2-complement is it really necessary to
       to set the entire path with a guard value */
    if (mo->model_type & GHMM_kSilentStates) {
        for (i = 0; i < len_path; i++) {
            state_seq[i] = -1;
        }
    }

    /* Precomputing the log(a_ij) and log(bj(ot)) */
    Viterbi_precompute(mo, o, len, v);

    /* Initialization, that is t = 0 */
    for (j = 0; j < mo->N; j++) {
        if (mo->s[j].pi == 0.0 || v->log_b[j][0] == +1) /* instead of 0, DBL_EPS.? */
            v->phi[j] = +1;
        else
            v->phi[j] = log(mo->s[j].pi) + v->log_b[j][0];
    }
    if (mo->model_type & GHMM_kSilentStates) {  /* could go into silent state at t=0 */
        __viterbi_silent(mo, t = 0, v);
    }

    /* t > 0 */
    for (t = 1; t < len; t++) {
        for (j = 0; j < mo->N; j++) {
            /* initialization of phi, psi */
            v->phi_new[j] = +1;
            v->psi[t][j] = -1;
        }

        for (k = 0; k < mo->N; k++) {
            /* Determine the maximum */
            /* max_phi = phi[i] + log_in_a[j][i] ... */
            if (!(mo->model_type & GHMM_kSilentStates) || !mo->silent[k]) {
                St = k;
                max_value = -DBL_MAX;
                v->psi[t][St] = -1;
                for (i = 0; i < mo->s[St].in_states; i++) {
                    if (v->phi[mo->s[St].in_id[i]] != +1 && v->log_in_a[St][i] != +1) {
                        value = v->phi[mo->s[St].in_id[i]] + v->log_in_a[St][i];
                        if (value > max_value) {
                            max_value = value;
                            v->psi[t][St] = mo->s[St].in_id[i];
                        }
                    }
                }

                /* No maximum found (that is, state never reached)
                   or the output O[t] = 0.0: */
                if (max_value == -DBL_MAX ||    /* and then also: (v->psi[t][j] == -1) */
                    v->log_b[St][t] == +1) {
                    v->phi_new[St] = +1;
                } else
                    v->phi_new[St] = max_value + v->log_b[St][t];

            }
        }                       /* complete time step for emitting states */

        /* First now replace the old phi with the new phi */
        for (j = 0; j < mo->N; j++) {
            v->phi[j] = v->phi_new[j];
            /*printf("\npsi[%d],%d, phi, %f\n", t, v->psi[t][j], v->phi[j]); */
        }

        /* complete time step for silent states */
        if (mo->model_type & GHMM_kSilentStates) {
            __viterbi_silent(mo, t, v);
        }
    }                           /* Next observation , increment time-step */

    /* Termination */
    /* for models with silent states we store the last state in the path at position 0.
       If there are no silent states we can use the correct position directly. */
    if (!(mo->model_type & GHMM_kSilentStates)) {
        state_seq_index = len_path - 1;
    } else {
        state_seq_index = 0;
    }

    max_value = -DBL_MAX;
    state_seq[state_seq_index] = -1;
    for (j = 0; j < mo->N; j++) {
        if (v->phi[j] != +1 && v->phi[j] > max_value) {
            max_value = v->phi[j];
            state_seq[state_seq_index] = j;
        }
    }
    if (max_value == -DBL_MAX) {
        /* Sequence can't be generated from the model! */
        *log_p = +1;
        if (!(mo->model_type & GHMM_kSilentStates)) {
            /* Backtracing doesn't work, insert -1 values in state_seq */
            for (t = len - 2; t >= 0; t--) {
                state_seq[t] = -1;
            }
        }
    } else {
        /* Backtracing, should put DEL path nicely */
        /* for models without silent states traceback is straightforward */
        if (!(mo->model_type & GHMM_kSilentStates)) {
            *log_p = max_value;
            lastemState = state_seq[len_path - 1];

            for (t = len - 2, i = len_path - 2; t >= 0; t--) {
                state_seq[i--] = v->psi[t + 1][lastemState];
                lastemState = v->psi[t + 1][lastemState];
            }
        }

        /* if there are silent states, we have to watch the length 
           of the viterbi path and realloc as needed */
        else {
            cur_len_path = 1;
            *log_p = max_value;
            lastemState = state_seq[0];

            for (t = len - 2, i = 1; t >= 0; t--) {
                /* if next state to be inserted into the path is silent we have to propagate up to the next emitting state */
                if (mo->silent[v->psi[t + 1][lastemState]]) {
                    St = v->psi[t + 1][lastemState];
                    /* fprintf(stderr, "t=%d:  DEL St=%d\n", t+1, St ); */

                    while (St != -1 && mo->silent[St]) {        /* trace-back up to the last emitting state */
                        /* fprintf(stderr, "***  t=%d:  DEL St=%d\n", t, St ); */
                        if (cur_len_path + 1 > len_path) {
                            /* we have to allocate more memory for state_seq. Memory is increased by the sequence length */
                            len_path = extend_int_array(state_seq, len_path, len);
                        }
                        state_seq[i++] = St;
                        St = v->psi[t][St];
                        cur_len_path++;
                    }

                    if (cur_len_path + 1 > len_path) {
                        /* we have to allocate more memory for state_seq. Memory is increased by the sequence length */
                        len_path = extend_int_array(state_seq, len_path, len);
                    }
                    state_seq[i++] = St;
                    lastemState = St;
                    cur_len_path++;

                } else {
                    if (cur_len_path + 1 > len_path) {
                        /* we have to allocate more memory for state_seq. Memory is increased by the sequence length */
                        len_path = extend_int_array(state_seq, len_path, len);
                    }

                    state_seq[i++] = v->psi[t + 1][lastemState];
                    lastemState = v->psi[t + 1][lastemState];
                    cur_len_path++;
                }
            }
        }
    }
    /* post-processing of the state path for models with silent states. 
       We have to realloc to the actual path length and reverse the order.
       The final element of the state path is marked with a -1 entry
       at the following array position */
    if (mo->model_type & GHMM_kSilentStates) {
        /* reallocating */
        if (cur_len_path + 1 != len_path) {
            ARRAY_REALLOC(state_seq, cur_len_path + 1);
            state_seq[cur_len_path] = -1;       /*end marker entry */
            len_path = cur_len_path + 1;
        }
        /* reversing order */
        for (i = 0; i < floor(cur_len_path / 2.0); i++) {
            k = state_seq[i];
            state_seq[i] = state_seq[cur_len_path - 1 - i];
            state_seq[cur_len_path - 1 - i] = k;
        }
    }

    /* Free the memory space */
    viterbi_free(&v, mo->N, len);
    return (state_seq);
  STOP:                        /* Label STOP from ARRAY_[CM]ALLOC */
    /* Free the memory space */
    viterbi_free(&v, mo->N, len);
    m_free(state_seq);
    return NULL;
#undef CUR_PROC
}                               /* ghmm_dmodel_viterbi */




/*============================================================================*/
double ghmm_dmodel_viterbi_logp(ghmm_dmodel * mo, int *o, int len, int *state_seq)
{
#define CUR_PROC "ghmm_dmodel_viterbi_logp"
    double log_p = 0.0;
    int *vpath = ghmm_dmodel_viterbi(mo, o, len, &log_p);
    m_free(vpath);
    return log_p;
#undef CUR_PROC
}                               /* ghmm_dmodel_viterbi_logp */
