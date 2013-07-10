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

/* based on  Speeding Up Bayesian HMM by the Four Russians Method
 *
 * M. Mahmud and A. Schliep
 *
 * In Algorithms in Bioinformatics, Springer Berlin / Heidelberg, 6833, 188–200, 2011. 
 * /

//=====================================================================
//====================precompute position==============================
//=====================================================================

//we need to be able to get the matrix M^B(X) where B is {0-M}^maxorder 
//and X is an obs sequence of length 1, 2, ..., R
//we precompute this and store index in stored pos
//currently B is only 0th order 

/* computes tupleSize and preposition
 * R: compression lenght
 * M: alphabet size
 * tupleSize: 0-M -> M, M + M^2, ...
 * preposition: preposition[i][j] = j*M^i */
void precomputeposition(int R, int M, int *tupleSize, int **preposition){
#define CUR_PROC "precomputeposition"
  //printf("m=%d\n", R);
  int i,j;
  int pw[R];
  tupleSize[0] = 0;
  tupleSize[1] = 0;
  pw[0] = 1;
  for(i = 1; i < R; i++){
    pw[i] = pw[i-1]*M;
    tupleSize[i+1] = tupleSize[i] + pw[i];
  }
  for(i = 0; i < R; i++){
    for(j = 0; j < M; j++){
      preposition[i][j] = j*pw[i];
      //printf("prepos %d %d = %d\n", i, j, preposition[i][j]);
    }
  }
#undef CUR_PROC
}

/* returns pos, the matrix correspoinding to the X = (o_start, o_start+1,..., obs_end) 
 * obs: observation sequence
 * start: start position in the sequence
 * end: end position in the sequence
 * tupleSize: 0-R -> M, M+M^2, ..
 * preposition: preposition[i][j] = j*M^i */
int position(int *obs, int start, int end, int *tupleSize, int **preposition){
#define CUR_PROC "position"
  int pos = tupleSize[end-start];
  int i = 0;
  for(; start<end; i++, start++){
    pos += preposition[i][obs[start]];
  }
  return pos;
#undef CUR_PROC
}

/* compute storedpos
 * R: compression length
 * T: length of obs 
 * storedpos: array to get matrix for each observation 0-limit -> M^B(X)
 * M: alphabet size
 * tupleSize: 0-M -> M, M + M^2, ...
 * preposition: preposition[i][j] = j*M^i */
void precomputeandstoreposition(int R, int T, int *obs, int *storedpos, int M, int *tupleSize, int **preposition){
#define CUR_PROC "precomputeandstoreposition"
    precomputeposition(R, M, tupleSize, preposition);
    int j, s, e, pos;
    storedpos[0] = position(obs, 0, R, tupleSize, preposition);
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
            pos = pos/M -1;//why not call position??
	    //printf("%d, e, pos = %d, position = %d\n",j,  pos,
		       	//position(obs, j, e, tupleSize, preposition));
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

void precomputeandstorepositionH(int R, int M, int order, int T, int *obs,
	       	int *storedpos,  int *storedfpos, int *tupleSize, int**preposition,
	       	int *tupleSizeH, int **prepositionH){
#define CUR_PROC "precomputeandstorepositionH"    
    int j, s, e, si ;
    int pos, fpos;

    precomputeposition(R, M, tupleSize, preposition);
    precomputeposition(order, M, tupleSizeH, prepositionH);

    storedpos[0] = position(obs, 0, R, tupleSize, preposition);
    storedfpos[0] = 0;

    si = 1 - order;
        if(si<0)
            si = 0;
    pos = position(obs, 1, R, tupleSize, preposition);
    fpos = position(obs, si, 1, tupleSizeH, prepositionH);
    storedpos[1] = pos;
    storedfpos[1] = fpos;
    for (j=2;j<R;j++){
        si = j - order;
        if(si<0)
            si = 0;
        pos = pos/M -1;
        fpos = position(obs, si, j, tupleSizeH, prepositionH);
        storedpos[j] = pos;
        storedfpos[j] = fpos;
    }

    s = R ;
    e = s + R ;
    si = s - order;
    if(si<0)
        si = 0;

    while (1){
        pos = position(obs, s, e, tupleSize, preposition);
        fpos = position(obs, si, s, tupleSizeH, prepositionH);
        storedpos[s] = pos ;
        storedfpos[s] = fpos;

        for (j=s+1;j<e;j++){
            si++;
            pos = pos/M -1;
            fpos = position(obs, si, j, tupleSizeH, prepositionH);
            storedpos[j] = pos;
            storedfpos[j] = fpos;
        }

        if (e == T)
            break;
        s += R;
        e += R;
        if (e>T)
            e = T;
        si = s- order ;
    }
#undef CUR_PROC
}

//=====================================================================
//=================precomputation matrix ==============================
//=====================================================================
//works todo dont alloc first M matrics of R (not used?)
//does all matrices some not needed?
/* R: compression length
 * mo: model
 * mats: M^B(X)
 * rmats: cdfs of mats used for sampling B3 */
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

//-----------------------HO-----------------------------------------

void recursivemats(int pos, int fpos, int a, int b, int R,
                int totalobs, int* obs, int **mflag,
                double ****mats, double *****rmats,
                int* storedpos, int *storedfpos, ghmm_dmodel *mo){

    if(a==b-1)
        return;

    if(mflag[pos][fpos])
        return;

    int npos = storedpos[a+1];//pos/AlphaSize - 1;
    int nfpos = storedfpos[a+1];//positionh(obs, fii, a+1);

    recursivemats(npos, nfpos, a+1, b, R, totalobs, obs, mflag, mats, rmats, storedpos, storedfpos, mo);
        
    int y = obs[a];
    double x;
    int j,k,l;
    double *p_rmats, *p_mats;
    for(j=0;j<mo->N;++j){
        for (k=0;k<mo->N;++k){
            p_rmats = &rmats[pos][fpos][j][k][0];
            p_mats = &mats[npos][nfpos][k][0];
            x = mats[y][fpos][0][j] * *p_mats ;
            *p_rmats = x ;
            p_rmats++;
            p_mats++;
            for (l=1; l<mo->N; ++l){
                x += mats[y][fpos][l][j] * *p_mats ;
                *p_rmats = x ;
                p_rmats++;
                p_mats++;
            }
            //printf("mats %d, %d, %d, %d = %f\n", pos, fpos, j, k, x);
            mats[pos][fpos][j][k] = x ;
        }
    }
    mflag[pos][fpos] = 1;
}


void lazyrecmats(int R, int totalobs, int *obs, int **mflag,
                double ****mats, double *****rmats, 
                int *storedpos, int *storedfpos, ghmm_dmodel* mo
                ){

    int s, e;
    for(s=0,e=R; ;){
        
        int pos = storedpos[s];//position(obs, s, e);
        int fpos = storedfpos[s];//positionh(obs, fi, s);

        recursivemats( pos, fpos, s, e, R, totalobs, obs, mflag, mats, rmats, storedpos, storedfpos, mo);

        s += R;
        e += R;
        if(s>=totalobs)
            break;
        if(e>totalobs)
            e = totalobs;
    }
}



void precomputedmatsH(int totalobs, int *obs, int R,
        double ****mats, double *****rmats,
        int **mflag,
        int *storedpos, int *storedfpos, ghmm_dmodel* mo){

    int i,j,k,l,m,n,e,pos,fpos;
    int dsize = (pow(mo->M, mo->maxorder+1)-1)/(mo->M-1);
    int size = (pow(mo->M, R+1)-1)/(mo->M-1);
    int read[(int)pow(mo->M, mo->maxorder)];
    int write[(int)pow(mo->M, mo->maxorder)];
    int *tmp1 = read;
    int *tmp2 = write;
    int *tmp3;
    for(i=0; i<size; i++) 
        for(j=0;j<dsize;j++)
            mflag[i][j] = 0;
    for(i=0;i<mo->M;i++)
        for(j=0;j<dsize;j++)
            mflag[i][j] = 1;

    
    fpos = 0;
    for(i=0; i <= mo->maxorder;i++)
    {
        n = 0;
        for(l = 0; l<pow(mo->M, i); l++)//fpos
        {
            for(pos = 0; pos < mo->M; pos++)
            { 
                
                for(j=0; j < mo->N; j++)
                {
                    for (k=0; k < mo->N; k++)
                    {
                        if(i == 0){
                            mo->emission_history = 0;
                            if(mo->order[k] ==0){
                                mats[pos][fpos][j][k] =
                                    ghmm_dmodel_get_transition(mo, j,k)*mo->s[k].b[pos];
                            }
                            else{
                                mats[pos][fpos][j][k] = 0;
                            } 
                        }
                        else{
                            mo->emission_history = tmp1[l];
                            e = get_emission_index(mo, k, pos, i);
                            mats[pos][fpos][j][k] =
                                ghmm_dmodel_get_transition(mo, j,k)*mo->s[k].b[e];
                        }
                        //printf("mats %d, %d, %d, %d, = %f\n", pos, fpos, j, k, 
                         //       mats[pos][fpos][j][k]);
                    }
                }
                update_emission_history_front(mo, m);
                tmp2[n] = mo->emission_history;                    
                n++;
            }
            fpos++;
        }
        tmp3 = tmp2;
        tmp2 = tmp1; 
        tmp1 = tmp3;
    }
    lazyrecmats(R, totalobs, obs, mflag, mats, rmats, storedpos, storedfpos, mo);
}


//===================================================================
//============= compressed sample state path ========================
//===================================================================
/* samples a dscrete distribution
 * seed: seed for uniform random sampler
 * distrubution: discrete distribution to sample
 * N: size of distrubution */
int samplebinsearch(int seed, double *distribution, int N){
#define CUR_PROC "samplebinsearch"
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
/* samples path Q_T~alpha_T
 * for m>=i<=2 and for s = (i-1)k, e = ik, Q_s ~ delta[s][q_e][*], 
 * for s<j<e, Q_j ~ R^o_j[Q_j-1][Q_e][*](o_j, o_j+1...e)
 * for j=1 to j=k-1 Q_j ~ M
 * seed: seed for random generator
 * T: length of observation sequence
 * obs: observeration sequence
 * fwds: forward variable
 * R: length of compression
 * mats: matrix form of foward algorithm used for compression
 * rmats: cdf of  mats used for sampling
 * states: states
 * storedpos: get matrix for observation
 * sneak: delta in pavels paper, cdfs of forwards
 * N: number of states */
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
            states[js] = samplebinsearch(seed, sneak[p+1][cs], N);
        }
        else{
            js = 0 ;
            pos = storedpos[1] ;
            states[js] = samplebinsearch(seed, sneak[p+1][cs], N);
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

void csamplestatepathH(int seed, int T, int *obs,
        double **fwds, int R, int N, 
        double ****mats, double *****rmats,
        int *states, int *storedpos, int *storedfpos,double ***sneak){

    double dist[N], total, rn;
    int j, s, e, md, p, l, r, m;
    int pos, cs, js, je, fpos, si;
    double *distribution;

    md = T%R;
    if (md == 0){
        md = R;
        p = T/R;
    }
    else{
        p = T/R + 1;
    }
    e = T;
    s = T - md;
    je = T - 2;
    
    states[e-1] = sample(seed, fwds[p], N);

    while ( s >= 0 ){
        p-- ;
        cs = states[e-1];
        if(s>0)
        {
            js = s - 1 ;
            pos = storedpos[s];
            fpos = storedfpos[s];
            states[js] = samplebinsearch(seed, sneak[p+1][cs], N);
        }
        else
        {
            js = 0;
            pos = storedpos[1];
            fpos = storedfpos[1];
            states[js] = samplebinsearch(seed, sneak[p+1][cs], N);
        }        

        for (;js<je;js++){
            distribution = rmats[pos][fpos][states[js]][cs];
            total = distribution[N-1];
            rn = ighmm_rand_uniform_cont(seed, total, 0);
            l = 0;
            r = N-1;
            while (1){
                m = (l+r)>>1;
                if (distribution[m] < rn)
                {
                    if(l==m){
                        states[js+1] = r;
                        break;
                    }
                    else
                        l = m;
                }
                else if (distribution[m] > rn)
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
            fpos = storedfpos[js+2];
        }
        s -= R;
        e = s + R;
        je = R + s - 2;
    }
}


//===================================================================
//==================compressed forwards==============================
//===================================================================
/* calulates forward variables and cdfs as a byproduct
 * totalobs: length of observation sequence
 * obs: observation sequence
 * R: compression length
 * fwds: forward variables
 * mats: matrix of forward varibles used for compression
 * storedpos: gets matrix for observation
 * sneak: cdf of mats columns */
void cforwards(int totalobs, int* obs, ghmm_dmodel *mo, int R, double **fwds, 
               double ***mats, int *storedpos, double ***sneak){
#define CUR_PROC "cforwards"
//why using matrix for forwards could use 2 vectors instead?     
    int i,j,k;
    double sum = 0, tv;
    int s, e ;
    int pos;
    //printf("first step\n");
    for (j=0;j<mo->N;j++){
        fwds[0][j] = mo->s[j].pi*mo->s[j].b[obs[0]];
        sum += fwds[0][j];
    }
    for (j=0;j<mo->N;j++){
        fwds[0][j] /= sum ;//XXX error checking
    }

    i = 1; 
    pos = storedpos[1];
    

    for (j=0;j<mo->N;j++){
        tv = fwds[0][0]*mats[pos][0][j];

        sneak[i][j][0] = tv;

        for (k=1;k<mo->N;k++){
            tv += fwds[0][k]*mats[pos][k][j];
            sneak[i][j][k] = tv;
        }
        fwds[i][j] = tv;
        sum += tv;
    }
    for(j=0;j<mo->N;j++)
        fwds[i][j] /= sum;

    i = 2;
    s = R;
    e = 2*R;
    
    while (1){

        pos = storedpos[s];
        sum = 0;
        for (j=0;j<mo->N;j++){
            tv = fwds[i-1][0]*mats[pos][0][j];
            sneak[i][j][0] = tv;
            for (k=1;k<mo->N;k++){
                tv += fwds[i-1][k]*mats[pos][k][j];
                sneak[i][j][k] = tv;
            }
            fwds[i][j] = tv;
            sum += tv;
        }

        for(j=0;j<mo->N;j++) 
            fwds[i][j] /= sum ;
        
        i++;
        s += R;
        if ( s >= totalobs )
            break;
    }        
#undef CUR_PROC
}

void cforwardsH(int T, int *obs, ghmm_dmodel *mo,
        int R, double **fwds,
        double ****mats, int *storedpos, int *storedfpos, double ***sneak
){

    int i,j,k;
    double sum = 0, tv;
    int s, e, si ;
    int pos, fpos;
    double *p_previous_fwds, *p_fwds;
    double *p_sneak;
    double *p_mats;

    //fist column of forwards
    for (j=0;j<mo->N;j++){
        if(mo->order[j] == 0)
            fwds[0][j] = mo->s[j].pi*mo->s[j].b[obs[0]];
        else
            fwds[0][j] = 0;
        sum += fwds[0][j];
    }
    for (j=0;j<mo->N;j++){
        fwds[0][j] /= sum ;
    }

    i = 1;
    
    
    p_sneak = &sneak[i][0][0];
    
    pos = storedpos[1];
    fpos = storedfpos[1];
    p_previous_fwds = fwds[i-1];
    p_fwds = fwds[i];
    sum = 0;
    for (j=0;j<mo->N;j++){        
        tv = p_previous_fwds[0]* mats[pos][fpos][0][j];
        sneak[i][j][0] = tv;
        for (k=1;k<mo->N;k++){
            tv += (p_previous_fwds[k]* mats[pos][fpos][k][j]) ;
            sneak[i][j][k] = tv;
        }
        p_fwds[j] = tv;
        sum += tv;
    }
    for(j=0;j<mo->N;j++)
         p_fwds[j] /= sum;

    i = 2 ;
    s = R ;
    e = s + R ;

    while (1){

        pos = storedpos[s];
        fpos = storedfpos[s];

        sum = 0;
        p_previous_fwds = p_fwds;
        p_fwds = fwds[i];
        
        for (j=0;j<mo->N;j++){            
            
            tv = p_previous_fwds[0]* mats[pos][fpos][0][j];
            sneak[i][j][0] = tv;

            for (k=1;k<mo->N;k++){
                tv += (p_previous_fwds[k]* mats[pos][fpos][k][j]) ;
                sneak[i][j][k] = tv;
            }

            p_fwds[j] = tv;
            sum += tv;

        }
        for(j=0;j<mo->N;j++){
            p_fwds[j] /= sum;
            //printf("fwd %d, %f\n",j, p_fwds[j]);
        }
        
        if (e == T)
            break;
        i++;
        s += R;
        e += R;
        if (e>T)
            e = T;
    }
}
//===================================================================
//===================== compressed gibbs ============================
//===================================================================
//no h.o, silence do error checking
/* runs the forward backward gibbs
 * mo: model
 * seed: seed 
 * obs: observation
 * totalobs: length of observation sequence
 * pA: prior for A
 * pB: prior for B
 * pPi: prior for pi
 * Q: states
 * R: length of compression */
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
//maybe alloc here instead and pass to step so no allocating every step
//should do precompute stuff here to not call it every iteration...doesnt change
/* runs the forward backward gibbs burnIn times
 * mo: model
 * seed: seed 
 * obs: observation
 * totalobs: length of observation sequence
 * pA: prior for A
 * pB: prior for B
 * pPi: prior for pi
 * Q: states
 * R: length of compression
 * burnIn: number of times to run forward backward gibbs */

void ghmm_dmodel_cfbgibbs(ghmm_dmodel* mo, int seed, int *obs, int totalobs, double **pA, double **pB, double *pPi, int *Q, int R, int burnIn){
    init_priors(&pA, &pB, &pPi, mo);
    for(;burnIn >= 0; burnIn--){
        if(mo->model_type & GHMM_kHigherOrderEmissions)
            ghmm_dmodel_cfbgibbstepH(mo, seed, obs, totalobs, pA, pB, pPi, Q, R);
        else 
            ghmm_dmodel_cfbgibbstep(mo, seed, obs, totalobs, pA, pB, pPi, Q, R);
    }
}


void ghmm_dmodel_cfbgibbstepH(ghmm_dmodel *mo, int seed, int *obs, int totalobs,
        double **pA, double **pB, double *pPi, int* Q, int R){
#define CUR_PROC "ghmm_dmodel_cfbgibbstepH"
        //mats, rmats, sneak, stored pos, fwds init here
        int shtsize = totalobs/R+2;
        double **fwds = ighmm_cmatrix_alloc(shtsize, mo->N);
        double ***sneak = ighmm_cmatrix_3d_alloc(shtsize, mo->N, mo->N);
        Q = malloc(sizeof(int)*totalobs);     
        
        //printf("precompute mat\n");
        double ****mats;
        double *****rmats;
        int i, j;
        int limit = (pow(mo->M, R+1)-1)/(mo->M-1);  
        int d = (pow(mo->M, mo->maxorder+1)-1)/(mo->M-1);
  	mats = malloc(sizeof(double*)*(limit+1));
  	rmats = malloc(sizeof(double*)*(limit+1));
        int **mflag = ighmm_dmatrix_alloc(limit, d);
        for(i = 0; i < limit+1; i++){
            mats[i] = ighmm_cmatrix_3d_alloc(d, mo->N, mo->N);
            rmats[i] = malloc(sizeof(double*)*d);
            for(j = 0; j < d; j++)
    	       rmats[i][j] = ighmm_cmatrix_3d_alloc(mo->N, mo->N, mo->N);
        }
        //printf("\nprecompute pos\n");
        int *tupleSize = malloc(sizeof(int)*(R+1));
        int *tupleSizeH = malloc(sizeof(int)*(mo->maxorder+1));
        int **preposition = ighmm_dmatrix_alloc(R, mo->M);
        int **prepositionH = ighmm_dmatrix_alloc(mo->maxorder, mo->M);
        int *storedpos = malloc(sizeof(int)*(totalobs+1));
        int *storedfpos = malloc(sizeof(int)*(totalobs+1));
        precomputeandstorepositionH(R, mo->M, mo->maxorder, totalobs, obs, storedpos, storedfpos, tupleSize, preposition, tupleSizeH, prepositionH);
  	precomputedmatsH(totalobs, obs, R, mats, rmats, mflag, storedpos, storedfpos, mo);
        //printf("\nforwards\n");
        cforwardsH(totalobs, obs, mo, R, fwds, mats, storedpos, storedfpos, sneak);
       
       // printf("\nsample\n");
        csamplestatepathH(seed, totalobs, obs, fwds, R, mo->N,  mats, rmats, Q, storedpos, storedfpos, sneak);

        //printf("\nupdate\n");
        updateH(seed, mo, totalobs, Q, obs, pA, pB, pPi);
        
       // printf("\nclean up\n");
        //clean up
        free(tupleSize);
        free(tupleSizeH);
        free(storedpos);
        free(storedfpos);
        ighmm_cmatrix_3d_free(&sneak, shtsize, mo->N);
        ighmm_dmatrix_free(&preposition, R);
        ighmm_dmatrix_free(&prepositionH, mo->maxorder);
        ighmm_cmatrix_free(&fwds, shtsize);
        for( i = 0 ; i < limit+1; i++){
           ighmm_cmatrix_3d_free(&mats[i],d, mo->N);
           for(j =0; j < d; j++)
              ighmm_cmatrix_3d_free(&(rmats[i][j]), mo->N, mo->N);
        } 
        free(rmats);
#undef CUR_PROC        
}
