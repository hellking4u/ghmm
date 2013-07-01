#include "cfbgibbs.h"
#include "matrix.h"
#include "ghmm.h"
#include "mes.h"
#include "ghmm_internals.h"
#include "sequence.h"
#include "rng.h"
#include "randvar.h"
#include <math.h>
#include "fbgibbs.h"

//=====================================================================
//====================precompute position==============================
//=====================================================================

//we map 0-limit to X, tupleSize allows us to jump to where the mapping to M^k starts
void precomputeposition(int R, int M, int *tupleSize, int **preposition){
#define CUR_PROC "precomputeposition"
  int i,j;
  int pw[R];
  tupleSize[1] = 0;
  pw[0] = 1;
  for(i = 1; i < R; i++){
    pw[i] = pw[i-1]*M;
    tupleSize[i+1] = tupleSize[i] + pw[i];
  }
  for(i = 0; i < R; i++){
    for(j = 0; j < M; j++){
      preposition[i][j] = j*pw[i];
    }
  }
#undef CUR_PROC
}

int position(int *obs, int start, int end, int *tupleSize, int **preposition){
#define CUR_PROC "position"
  int pos = tupleSize[end-start];
  int i = 0;
  for(; start<end; i++< start++){
    pos += preposition[i][obs[start]];
  }
  return pos;
#undef CUR_PROC
}

void precomputeandstoreposition(int R, int T, int *obs, int *storedpos, int M, int *tupleSize, int **preposition){
#define CUR_PROC "precomputeandstoreposition"
    precomputeposition(R, M, tupleSize, preposition);
    int j, s, e, pos;
    pos = position(obs, 1, R, tupleSize, preposition);
    storedpos[1] = pos;    
    for(j=2;j<R;j++){
        pos = pos/M -1;
        storedpos[j] = pos;
    }
    
    s = R;
    e = s + R;    
    while (1){
        pos = position(obs, s, e, tupleSize, preposition);
        storedpos[s] = pos;
        for(j=s+1;j<e;j++){
            pos = pos/M -1;
            storedpos[j] = pos;
        }

        if (e == T)
            break;
        s += R;
        e += R;
        if (e>T)
            e = T;
    }
#undef CUR_PROC
}

//=====================================================================
//=================precomputation matrix ==============================
//=====================================================================
//works todo dont alloc first M matrics of R (not used)

//pavel four russians paper
//mo is the model
//mats is M
//rmats is R
void precompute(int R, ghmm_dmodel *mo, double ***mats, double**** rmats){
#define CUR_PROC "precompute"
  int i, j, k;
  int limit = (pow(mo->M, R+1)-1)/(mo->M-1) -1;  
  //mats[i][j][k] i = obs; j,k indice of matrix   M(obs)_ik 
  for(i = 0; i < mo->N; i++){
    for(j = 0; j < mo->N; j++){
      for(k = 0; k < mo->M; k++){
        mats[k][i][j] = ghmm_dmodel_get_transition(mo, i,j)*mo->s[j].b[k];
        //printf("mats(%d, %d, %d) = %f \n",i,j,k, mats[i][j][k]);
      }
    }  
  }
  
 int pos = mo->M;
 int fpos = 0;
 int y = 0;
 double x;
 while(pos < limit){ 
     for(j = 0; j < mo->N; j++){
         for(k = 0; k < mo->N; k++){
             x = mats[y][j][0] * mats[fpos][0][k];
             rmats[pos][j][k][0] = x;
             //printf("rmats(%d, %d, %d, %d) = %f \n",pos,j,k,0, rmats[pos][j][k][0]);
             for(i = 1; i < mo->N; i++){
                 x += mats[y][j][i] * mats[fpos][i][k]; 
                 rmats[pos][j][k][i] = x;
                 //printf("rmats(%d, %d, %d, %d) = %f \n",pos,j,k,i, rmats[pos][j][k][i]);
             }
                 mats[pos][j][k] = x;
	//printf("mats(%d, %d, %d) = %f \n",pos,j,k, mats[pos][j][k]);
         }
    }
    pos++;
    y++;
    if(y == mo->M){
        y = 0;
        fpos = pos/mo->M - 1; 
        //printf("fpos = %d\n", fpos); 
    } 
  }
#undef CUR_PROC
}
//===================================================================
//============= compressed sample state path ========================
//===================================================================
int samplestatebinsearch(int seed, double *distribution, int N){
#define CUR_PROC "samplestatebinsearch"
    double total = distribution[N-1];
    double rn = ighmm_rand_uniform_cont(seed, total, 0);
    int l = 0;
    int r = N-1;
    int m;
    while (1){
        m = (l+r)>>1;
        if (distribution[m] < rn)
        {
            if(l==m)
                return r;
            else
                l = m;
        }
        else if (distribution[m] > rn)
        {
            if(l==m)
                return l;
            else
                r = m;
        }
        else
        {
           return m;
        }
    }
#undef CUR_PROC
}
void csamplestatepath(int seed, int T, int *obs,
        double **fwds, int R,
        double ***mats, double ****rmats,
        int *states, int* storedpos, double ***sneak, int N){
#define CUR_PROC "csamplestatepath"
    double dist[N], total, rn, tv;
    double *distribution;
    int pos, cs, js, je;
    int p, s, e, l, r, m ;
    int md = T%R ;
    
    if (md == 0){
        md = R;
        p = T/R;
    }
    else{
        p = T/R + 1;
    }
    e = T;
    s = T - md;

    states[e-1] = sample(seed, fwds[p], N);
           
    while (s>=0){
        p--;
        cs = states[e-1];

        if(s>0){
            js = s - 1;
            pos = storedpos[s];                       
            states[js] = samplestatebinsearch(seed, sneak[p+1][cs], N);
        }
        else{
            js = 0 ;
            pos = storedpos[1] ;
            states[js] = samplestatebinsearch(seed, sneak[p+1][cs], N);
        }
        
        je = md + s -2;
        
        for (;js<je;js++){
            distribution = rmats[pos][states[js]][cs];
            total = distribution[N-1];
            rn = ighmm_rand_uniform_cont(seed, total, 0);

            l = 0;
            r = N-1;
            while (1){
                m = (l+r)>>1;
                tv = distribution[m];
                if ( tv < rn)
                {
                    if(l==m){
                        states[js+1] = r;
                        break;
                    }
                    else
                        l = m;
                }
                else if ( tv > rn)
                {
                    if(l==m){
                        states[js+1] =l;
                        break;
                    }
                    else
                        r = m;
                }
                else{
                    states[js+1] = m;
                    break;
                }
            }

            pos = storedpos[js+2];

        }

        md = R;
        s -= md;
        e = s + md;
    }
#undef CUR_PROC
}
//===================================================================
//==================compressed forwards==============================
//===================================================================
void cforwards(int totalobs, int* obs, ghmm_dmodel *mo, long R, double **fwds, 
               double ***mats, int *storedpos, double ***sneak){
#define CUR_PROC "cforwards"
//why using matrix for forwards could use 2 vectors instead?     
    int i,j,k;
    double sm = 0, tv;
    int s, e ;
    int prtobs;
    //printf("first step\n");
    for (j=0;j<mo->N;j++){
        fwds[0][j] = mo->s[j].pi*mo->s[j].b[obs[0]];
        sm += fwds[0][j];
    }
    for (j=0;j<mo->N;j++){
        fwds[0][j] /= sm ;//XXX error checking
    }

    i = 1; 
    prtobs = storedpos[1];
    

    for (j=0;j<mo->N;j++){
        tv = fwds[0][0]*mats[prtobs][0][j];

        sneak[i][j][0] = tv;

        for (k=1;k<mo->N;k++){
            tv += fwds[0][k]*mats[prtobs][k][j];
            sneak[i][j][k] = tv;
        }
        fwds[i][j] = tv;
        sm += tv;
    }
    for(j=0;j<mo->N;j++)
        fwds[i][j] /= sm;

    i = 2;
    s = R;
    e = 2*R;
    
    while (1){

        prtobs = storedpos[s];
        sm = 0;
        for (j=0;j<mo->N;j++){
            tv = fwds[i-1][0]*mats[prtobs][0][j];
            sneak[i][j][0] = tv;
            for (k=1;k<mo->N;k++){
                tv += fwds[i-1][k]*mats[prtobs][k][j];
                sneak[i][j][k] = tv;
            }
            fwds[i][j] = tv;
            sm += tv;
        }

        for(j=0;j<mo->N;j++) 
            fwds[i][j] /= sm ;
        
        i++;
        s += R;
        if ( s >= totalobs )
            break;
    }        
#undef CUR_PROC
}

//===================================================================
//===================== compressed gibbs ============================
//===================================================================
//no h.o, silence do error checking
void ghmm_dmodel_cfbgibbstep(ghmm_dmodel *mo, int seed, int *obs, int totalobs,
        double **pA, double **pB, double *pPi, int* Q, int R){
#define CUR_PROC "ghmm_dmodel_cfbgibbstep"
        //mats, rmats, sneak, stored pos, fwds init here

        int shtsize = totalobs/R+2;
        double **fwds = ighmm_cmatrix_alloc(shtsize, mo->N);
        double ***sneak = ighmm_cmatrix_3d_alloc(shtsize, mo->N, mo->N);
        
        
        //printf("precompute mat\n");
        double ***mats;
        double ****rmats;
        int i;
        int limit = (pow(mo->M, R+1)-1)/(mo->M-1) -1;  
  	mats = ighmm_cmatrix_3d_alloc(limit, mo->N, mo->N);
  	rmats = malloc(sizeof(double*)*limit);
  	for(i = 0; i < limit; i++)
    	    rmats[i] = ighmm_cmatrix_3d_alloc(mo->N, mo->N, mo->N);
        precompute(R, mo, mats, rmats);

        
        //printf("\nprecompute pos\n");
        int *tupleSize = malloc(sizeof(int)*(R+1));
        int **preposition = ighmm_dmatrix_alloc(R, mo->M);
        int *storedpos = malloc(sizeof(int)*(totalobs+1));
        precomputeandstoreposition(R, totalobs, obs, storedpos, mo->M, tupleSize, preposition);

        
        //printf("\nforwards\n");
        cforwards(totalobs, obs, mo, R, fwds, mats, storedpos, sneak);
        
       // printf("\nsample\n");
        csamplestatepath(seed, totalobs, obs, fwds, R, mats, rmats, Q, storedpos, sneak, mo->N);
        
        //printf("\nupdate\n");
        update(seed, mo, totalobs, Q, obs, pA, pB, pPi);
        
       // printf("\nclean up\n");
        //clean up
        free(tupleSize);
        free(storedpos);
        ighmm_cmatrix_3d_free(&sneak, shtsize, mo->N);
        ighmm_dmatrix_free(&preposition, R);
        ighmm_cmatrix_free(&fwds, shtsize);
        ighmm_cmatrix_3d_free(&mats,limit, mo->N);
        for( i = 0 ; i < limit; i++)
           ighmm_cmatrix_3d_free(&(rmats[i]), mo->N, mo->N);
        free(rmats);
#undef CUR_PROC        
}
void ghmm_dmodel_cfbgibbs(ghmm_dmodel* mo, int seed, int *obs, int totalobs, double **pA, double **pB, double *pPi, int *Q, int R, int burnIn){
	for(;burnIn >= 0; burnIn--){
            ghmm_dmodel_cfbgibbstep(mo, seed, obs, totalobs, pA, pB, pPi, Q, R);
        }

}
