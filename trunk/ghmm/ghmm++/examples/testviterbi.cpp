
#include <ghmm/matrix.h>
#include <ghmm/vector.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/sdmodel.h>
#include "ghmm/mes.h"


#include <stdio.h>
#include <stdlib.h>
#include "ghmm++/GHMM.h"
#include "ghmm++/GHMM_Computation.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

static double **mat_d_alloc(int n, int m)
{
  int i,j;
  double** matout = new double*[n];
  for (i=0;i<n;i++)
    {
      matout[i] = new double[m];
    }
  return matout;
}

static void mat_d_free(double ***mat, int n)
{
  int i;
  for (i=0;i<n;i++)
    {
      delete [] (*mat)[i];
    }
  delete *mat;
}

int seq_d_class(const double *O, int index, double *osum) {
  printf(" c:%i ",index);
  if (index>5)
    return 1;
  else
    return 0; 
}

int* reduce_vec(int** vec,int n)
{
  int i = 0;
  int k = 0;
  int* vecout;
  while (((*vec)[k]!=-1) && (k<n))
      k++;
  vecout = new int[k];
  for (i=0;i<k;i++)
    vecout[i] = (*vec)[i];
  delete *vec;
  return vecout;
}

void create_model(sdmodel* sdm,double*** matrans,double** maems,double* vecpi,int n,int m,int cos)
{
  int i,j,k,l,c;
  double** out_probs;
  double** in_probs;
  int* out_ids;
  int* in_ids;

  assert(sdm);

  sdm->get_class = seq_d_class;
  sdm->N = n; //number of states
  sdm->M = m; //length of alphabet
  sdm->cos = cos; //number of classes
  sdm->prior = -1;
  sdm->s = new sdstate[n];
  assert(sdm->s);

  matrix_d_print(stdout,maems,n,m,"\t"," ","\n");
  for (i=0;i<cos;i++)
    {
      matrix_d_print(stdout,matrans[i],n,m,"\t"," ","\n");
    }
  
  for (i=0;i<n;i++)
    {
      out_ids = new int[n];
      in_ids  = new int[n];
      out_ids[k=0] = -1;
      in_ids[l=0] = -1;

      for (j=0;j<n;j++)
	{
	  for (c=0;c<cos;c++)
	    {
	      //assigning out_states
	      if ((out_ids[k]==-1) && (matrans[c][i][j]>0))
		  out_ids[k] = j;
	      //assigning in_states
	      if ((in_ids[l]==-1) && (matrans[c][j][i]>0))
		  in_ids[l] = j;
	    }
	  if (out_ids[k]!=-1)
	    {
	      k++;
	      if (k<n)
		out_ids[k] = -1;
	    }
	  if (in_ids[l]!=-1)
	    {
	      l++;
	      if (l<n)
		in_ids[l] = -1;
	    }
	}

      out_ids = reduce_vec(&out_ids,n);
      in_ids = reduce_vec(&in_ids,n);
      sdm->s[i].pi = vecpi[i];
      sdm->s[i].b = maems[i];
      sdm->s[i].out_states = k;
      sdm->s[i].in_states = l;

      out_probs = mat_d_alloc(cos,k);
      in_probs  = mat_d_alloc(cos,l);

      //out_states
      for (j=0;j<k;j++)
	  for(c=0;c<cos;c++)
	      out_probs[c][j] = matrans[c][i][out_ids[j]];
      //in_states
      for (j=0;j<l;j++)
	  for(c=0;c<cos;c++)
	      in_probs[c][j] = matrans[c][out_ids[j]][i];

      sdm->s[i].out_id = out_ids;
      sdm->s[i].in_id = in_ids;
      sdm->s[i].out_a = out_probs;
      sdm->s[i].in_a = in_probs;
      
    }
}

double** to_2dyn(double* mat,int n,int m)
{
  int i,j;
  double** matout = mat_d_alloc(n,m);
  for (i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
	{
	  matout[i][j] = mat[i*m+j];
	}
    }
  //matrix_d_print(stdout,matout,n,m,"\t"," ","\n");
  return matout;
}

double*** to_3dyn(double* mat,int n,int m, int k)
{
  int i,j,l;
  double*** matout = new double**[k];
  for (i=0;i<k;i++)
    {
      matout[i] = to_2dyn(&mat[i*n*m],n,m);
      //matrix_d_print(stdout,matout[i],n,m,"\t"," ","\n");
    }
  return matout;
}


void mat3_free(double ****matrix, long levels, long rows)
{
  int i;
  for (i=levels-1;i>=0;i--)
    {
      mat_d_free((double***)&((*matrix)[i]),rows);
    }
  delete *matrix;
}

int main(int argc,char* argv)
{
  int n = 2; //number of states
  int m = 4; //length of alphabet
  int cos = 2; //number of output classes
  double mt[2][2][2] = {{{1.0, 0.0},{0.0,1.0}},{{0.0,1.0},{0.0,1.0}}};
  double*** matrans = to_3dyn((double*)mt,n,n,cos);
  double me[2][4] = {{0.4,0.3,0.3,0.0},{0.0,0.0,0.0,1.0}};
  double** maems = to_2dyn((double*)me,n,m);
  double vecpi[2] = {1.0,0.0};
  sequence_t* seqs;

  Viterbi* vbiobj;

  gsl_rng_init(); //Init - important!!

  sdmodel* sdm = new sdmodel;
  create_model(sdm,matrans,maems,vecpi,n,m,cos);


  printf("Klasse 1:\n");
  sdmodel_Ak_print( stdout, sdm ,0,"\t"," ","\n");
  printf("Klasse 2:\n");
  sdmodel_Ak_print( stdout, sdm ,1,"\t"," ","\n");
  sdmodel_B_print(stdout,sdm,"\t"," ","\n");

  seqs = sdmodel_generate_sequences(sdm,1,10,10,10);  
  sequence_print(stdout,seqs);

  vbiobj = new Viterbi(sdm, 10);
  double logp = vbiobj->Viterbi_runme(sdm, seqs->seq[9], 10);

  printf("\n\n Viterbi path %f", logp );
  vbiobj->print_path(10, "\t");

  mat3_free(&matrans,cos,n);
  matrix_d_free(&maems,n);
  sdmodel_free(&sdm);

  delete vbiobj;

}
