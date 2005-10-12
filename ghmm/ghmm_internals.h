#ifndef GHMM_INTERNALS_H
#define GHMM_INTERNALS_H

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/


/*==============  Memory allocation macros for mes-functions  ===============*/
#ifndef ARRAY_MALLOC
#define ARRAY_MALLOC(ptr, entries) { if (!((ptr) = ighmm_malloc(sizeof(*(ptr)) * (entries)))) \
                                         {mes_proc(); goto STOP;}                           \
				   }
#endif

#ifndef ARRAY_CALLOC
#define ARRAY_CALLOC(ptr, entries) { if (!((ptr) = ighmm_calloc(sizeof(*(ptr)) * (entries))))   \
                                         {mes_proc(); goto STOP;}                             \
				   }
#endif

#ifndef ARRAY_REALLOC
#define ARRAY_REALLOC(ptr, entries) { if (ighmm_realloc ((void**)&(ptr), sizeof(*(ptr))*(entries))) \
                                          {mes_proc(); goto STOP;}                                \
				    }
#endif


/*==============  declarations for root_finder.c  ===========================*/
/**
   brent root finding algorithm.
   wrapps this functioncall to the gsl 1D-root finding interface
   @author Achim G\"adke
   @param x1 lower bracket value
   @param x2 upper bracket value
   @param tol tolerance for iteration
   @param A 1st extra parameter
   @param B 2nd extra parameter
   @param eps 3rd extra parameter
   @return root
 */
  double ghmm_zbrent_AB (double (*func) (double, double, double, double),
                    double x1, double x2, double tol, double A, double B,
                    double eps);


/*==============  declarations for gauss_tail.c  ============================*/
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


/*==============  numeric functions  ========================================*/
/**
   Calculates the logarithm of the sum of the probabilities whose logarithms are
   stored in the given array
   @return log of sum of exp(a[i])
   @param a:          array of logarithms of probabilities (a[i] < 0 for all i)
   @param length:     length of a
*/
  double ighmm_cvector_log_sum (double *a, int N);


#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif
