#include "nr.h"
#include <math.h>
#include <stdlib.h>

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b), (maxarg1) > (maxarg2) ?\
		   (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


#define NR_END 1
#define FREE_ARG char*

static void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);  
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}



/* EIGENBAU: */
static void nrwarning(char error_text[])
{
  fprintf(stderr,"Numerical Recipes warning...\n");
  fprintf(stderr,"%s\n",error_text);
}


/*============================================================================*/
#define ITMAX 100
#define EPS 3.0e-8

/* Eigene Aenderungen: func um 3 Parameter A, B, eps erweitert, die auch mit 
   uebergeben werden (vgl. Anwendung in creestimate bzw. sreestimate);
   (BK, 2.9.99) 
   aus float jeweils double gemacht.
   (BK, 12.10.99)
*/
double zbrent_AB(double (*func)(double, double, double, double),
		 double x1, double x2, double tol, double A, double B, 
		 double eps)
{
  int iter;
  double a=x1,b=x2,c=x2,d,e,min1,min2;
  double fa=(*func)(a, A, B, eps), fb=(*func)(b, A, B, eps),
    fc,p,q,r,s,tol1,xm;
  

  /* TEST:
  int i=0;
  while (fa < 0.0 && a < b && i<5) { 
    a *= 0.9999; fa=(*func)(a, A, B, eps); i++;
  }*/
  if (fa < 0.0 && fb < 0.0) {
    //printf("A = %e, B = %e, mue_l = %e, mue_r = %e\n", A, B, a, b);
    //printf("pmue(mue_left) < 0.0! -> mue := mue_left\n");
    //nrwarning("pmue(mue_l) < 0.0!");
    return a;
  }
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    nrerror("Root must be bracketed in zbrent");
  }
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) {
      //printf("#iter=%d, tol1=%.8f\n",iter,tol1);
      return b;
    }
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=(*func)(b, A, B, eps);
  }
  /* never get here */
  fprintf(stderr, "min1 %e, min2 %e\n", min1, min2);
  fprintf(stderr,"x1 %e, x2 %e, tol %e, A %e, B%e\n", x1, x2, tol, A, B);
  fprintf(stderr,"a %e, b %e, c %e, d %e, e %e, fa %e, fb %e, fc %e\n", 
	  a, b, c, d, e, fa, fb, fc);
  nrerror("Maximum number of iterations exceeded in zbrent");
  return 0.0;
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
#define MAXIT 100

double rtsafe(void (*funcd)(double, double *, double *), double x1, double x2,
	      double xacc)
{
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;
  
  (*funcd)(x1,&fl,&df);
  (*funcd)(x2,&fh,&df);
  /* if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
       nrerror("Root must be bracketed in rtsafe"); */

  /* Modifikation: Falls keine Nullstelle: kleinster Rand */
  if (fl > 0.0 && fh > 0.0) {
    /* nrwarning("both interval function values > 0!"); */
    if (fl < fh) return x1;
    else return x2;
  }
  else if (fl < 0.0 && fh < 0.0) {
    /* nrwarning("both interval function values < 0!"); */
    if (fl > fh) return x1;
    else return x2;
  }

  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts,&f,&df);
  for (j=1;j<=MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
	|| (fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts;
    (*funcd)(rts,&f,&df);
    if (f < 0.0)
      xl=rts;
    else
      xh=rts;
  }
  nrerror("Maximum number of iterations exceeded in rtsafe");
  return -1.0;
}
#undef MAXIT
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */



/*============================================================================*/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void lnsrch(int n, double xold[], double fold, double g[], double p[], 
	    double x[], double *f, double stpmax, int *check, 
	    double (*func)(double []))
{
  int i;
  double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
    test,tmplam;
  
  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for (i=1;i<=n;i++) {
    temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*func)(x);
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
      return;
    } else if (*f <= fold+ALF*alam*slope) return;
    else {
      if (alam == 1.0)
	tmplam = -slope/(2.0*(*f-fold-slope));
      else {
	rhs1 = *f-fold-alam*slope;
	rhs2=f2-fold2-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc<0.0) nrerror("Roundoff problem in lnsrch.");
	  else tmplam=(-b+sqrt(disc))/(3.0*a);
	}
	if (tmplam>0.5*alam)
	  tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=FMAX(tmplam,0.1*alam);
  }
}
#undef ALF
#undef TOLX
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */




/*============================================================================*/

#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n); \
free_matrix(hessin,1,n,1,n);free_vector(hdg,1,n);free_vector(g,1,n); \
free_vector(dg,1,n);

void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double []), void (*dfunc)(double [], double []))
{
	void lnsrch(int n, double xold[], double fold, double g[], double p[], 
		    double x[], double *f, double stpmax, int *check, 
		    double (*func)(double []));
	int check,i,its,j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
	double *dg,*g,*hdg,**hessin,*pnew,*xi;

	dg=vector(1,n);
	g=vector(1,n);
	hdg=vector(1,n);
	hessin=matrix(1,n,1,n);
	pnew=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,g);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
		fp = *fret;
		for (i=1;i<=n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i];
		(*dfunc)(p,g);
		test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
			FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac*fac > EPS*sumdg*sumxi) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
				}
			}
		}
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	nrerror("too many iterations in dfpmin");
	FREEALL
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */




/*============================================================================*/

#define EPS 1.0e-4

void fdjac(int n, double x[], double fvec[], double **df,
	void (*vecfunc)(int, double [], double []))
{
	int i,j;
	double h,temp,*f;

	f=vector(1,n);
	for (j=1;j<=n;j++) {
		temp=x[j];
		h=EPS*fabs(temp);
		if (h == 0.0) h=EPS;
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(n,x,f);
		x[j]=temp;
		for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
	}
	free_vector(f,1,n);
}
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/

extern int nn;
extern double *fvec;
extern void (*nrfuncv)(int n, double v[], double f[]);

double fmin(double x[])
{
	int i;
	double sum;

	(*nrfuncv)(nn,x,fvec);
	for (sum=0.0,i=1;i<=nn;i++) sum += SQR(fvec[i]);
	return 0.5*sum;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/

void rsolv(double **a, int n, double d[], double b[])
{
	int i,j;
	double sum;

	b[n] /= d[n];
	for (i=n-1;i>=1;i--) {
		for (sum=0.0,j=i+1;j<=n;j++) sum += a[i][j]*b[j];
		b[i]=(b[i]-sum)/d[i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/

void qrdcmp(double **a, int n, double *c, double *d, int *sing)
{
  int i,j,k;
  double scale=0.0,sigma,sum,tau;
  
  *sing=0;
  for (k=1;k<n;k++) {
    for (i=k;i<=n;i++) scale=FMAX(scale,fabs(a[i][k]));
    if (scale == 0.0) {
      *sing=1;
      c[k]=d[k]=0.0;
    } else {
      for (i=k;i<=n;i++) a[i][k] /= scale;
      for (sum=0.0,i=k;i<=n;i++) sum += SQR(a[i][k]);
      sigma=SIGN(sqrt(sum),a[k][k]);
      a[k][k] += sigma;
      c[k]=sigma*a[k][k];
      d[k] = -scale*sigma;
      for (j=k+1;j<=n;j++) {
	for (sum=0.0,i=k;i<=n;i++) sum += a[i][k]*a[i][j];
	tau=sum/c[k];
	for (i=k;i<=n;i++) a[i][j] -= tau*a[i][k];
      }
    }
  }
  d[n]=a[n][n];
  if (d[n] == 0.0) *sing=1;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/

void rotate(double **r, double **qt, int n, int i, double a, double b)
{
	int j;
	double c,fact,s,w,y;

	if (a == 0.0) {
		c=0.0;
		s=(b >= 0.0 ? 1.0 : -1.0);
	} else if (fabs(a) > fabs(b)) {
		fact=b/a;
		c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
		s=fact*c;
	} else {
		fact=a/b;
		s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
		c=fact*s;
	}
	for (j=i;j<=n;j++) {
		y=r[i][j];
		w=r[i+1][j];
		r[i][j]=c*y-s*w;
		r[i+1][j]=s*y+c*w;
	}
	for (j=1;j<=n;j++) {
		y=qt[i][j];
		w=qt[i+1][j];
		qt[i][j]=c*y-s*w;
		qt[i+1][j]=s*y+c*w;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/

void qrupdt(double **r, double **qt, int n, double u[], double v[])
{
	void rotate(double **r, double **qt, int n, int i, double a, double b);
	int i,j,k;

	for (k=n;k>=1;k--) {
		if (u[k]) break;
	}
	if (k < 1) k=1;
	for (i=k-1;i>=1;i--) {
		rotate(r,qt,n,i,u[i],-u[i+1]);
		if (u[i] == 0.0) u[i]=fabs(u[i+1]);
		else if (fabs(u[i]) > fabs(u[i+1]))
			u[i]=fabs(u[i])*sqrt(1.0+SQR(u[i+1]/u[i]));
		else u[i]=fabs(u[i+1])*sqrt(1.0+SQR(u[i]/u[i+1]));
	}
	for (j=1;j<=n;j++) r[1][j] += u[1]*v[j];
	for (i=1;i<k;i++)
		rotate(r,qt,n,i,r[i][i],-r[i+1][i]);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/

#define MAXITS 200
#define EPS 1.0e-7
#define TOLF 1.0e-4
#define TOLX EPS
#define STPMX 100.0
#define TOLMIN 1.0e-6
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
	free_vector(w,1,n);free_vector(t,1,n);free_vector(s,1,n);\
	free_matrix(r,1,n,1,n);free_matrix(qt,1,n,1,n);free_vector(p,1,n);\
	free_vector(g,1,n);free_vector(fvcold,1,n);free_vector(d,1,n);\
	free_vector(c,1,n);return;}

int nn;
double *fvec;
void (*nrfuncv)(int n, double v[], double f[]);

void broydn(double x[], int n, int *check,
	    void (*vecfunc)(int, double [], double [])) {
  void fdjac(int n, double x[], double fvec[], double **df,
	     void (*vecfunc)(int, double [], double []));
  double fmin(double x[]);
  void lnsrch(int n, double xold[], double fold, double g[], double p[], 
	      double x[], double *f, double stpmax, int *check, 
	      double (*func)(double []));
  void qrdcmp(double **a, int n, double *c, double *d, int *sing);
  void qrupdt(double **r, double **qt, int n, double u[], double v[]);
  void rsolv(double **a, int n, double d[], double b[]);
  int i,its,j,k,restrt,sing,skip;
  double den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
  double *g,*p,**qt,**r,*s,*t,*w,*xold;
  
  c=vector(1,n);
  d=vector(1,n);
  fvcold=vector(1,n);
  g=vector(1,n);
  p=vector(1,n);
  qt=matrix(1,n,1,n);
  r=matrix(1,n,1,n);
  s=vector(1,n);
  t=vector(1,n);
  w=vector(1,n);
  xold=vector(1,n);
  fvec=vector(1,n);
  nn=n;
  nrfuncv=vecfunc;
  f=fmin(x);
  test=0.0;
  for (i=1;i<=n;i++)
    if (fabs(fvec[i]) > test)test=fabs(fvec[i]);
  if (test<0.01*TOLF) FREERETURN
			for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  restrt=1;
  for (its=1;its<=MAXITS;its++) {
    if (restrt) {
      fdjac(n,x,fvec,r,vecfunc);
      qrdcmp(r,n,c,d,&sing);
      if (sing) nrerror("singular Jacobian in broydn");
      for (i=1;i<=n;i++) {
	for (j=1;j<=n;j++) qt[i][j]=0.0;
	qt[i][i]=1.0;
      }
      for (k=1;k<n;k++) {
	if (c[k]) {
	  for (j=1;j<=n;j++) {
	    sum=0.0;
	    for (i=k;i<=n;i++)
	      sum += r[i][k]*qt[i][j];
	    sum /= c[k];
	    for (i=k;i<=n;i++)
	      qt[i][j] -= sum*r[i][k];
	  }
	}
      }
      for (i=1;i<=n;i++) {
	r[i][i]=d[i];
	for (j=1;j<i;j++) r[i][j]=0.0;
      }
    } else {
      for (i=1;i<=n;i++) s[i]=x[i]-xold[i];
      for (i=1;i<=n;i++) {
	for (sum=0.0,j=i;j<=n;j++) sum += r[i][j]*s[j];
	t[i]=sum;
      }
      skip=1;
      for (i=1;i<=n;i++) {
	for (sum=0.0,j=1;j<=n;j++) sum += qt[j][i]*t[j];
	w[i]=fvec[i]-fvcold[i]-sum;
	if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i]))) skip=0;
	else w[i]=0.0;
      }
      if (!skip) {
	for (i=1;i<=n;i++) {
	  for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*w[j];
	  t[i]=sum;
	}
	for (den=0.0,i=1;i<=n;i++) den += SQR(s[i]);
	for (i=1;i<=n;i++) s[i] /= den;
	qrupdt(r,qt,n,t,s);
	for (i=1;i<=n;i++) {
	  if (r[i][i] == 0.0) nrerror("r singular in broydn");
	  d[i]=r[i][i];
	}
      }
    }
    for (i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
      g[i]=sum;
    }
    for (i=n;i>=1;i--) {
      for (sum=0.0,j=1;j<=i;j++) sum += r[j][i]*g[j];
      g[i]=sum;
    }
    for (i=1;i<=n;i++) {
      xold[i]=x[i];
      fvcold[i]=fvec[i];
    }
    fold=f;
    for (i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
      p[i] = -sum;
    }
    rsolv(r,n,d,p);
    lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
    test=0.0;
    for (i=1;i<=n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) {
      *check=0;
      FREERETURN
	}
    if (*check) {
      if (restrt) FREERETURN
		    else {
		      test=0.0;
		      den=FMAX(f,0.5*n);
		      for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
			if (temp > test) test=temp;
		      }
		      if (test < TOLMIN) FREERETURN
					   else restrt=1;
		    }
    } else {
      restrt=0;
      test=0.0;
      for (i=1;i<=n;i++) {
	temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
	if (temp > test) test=temp;
      }
      if (test < TOLX) FREERETURN
			 }
  }
  nrerror("MAXITS exceeded in broydn");
  FREERETURN
}
#undef MAXITS
#undef EPS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort(unsigned long n, double arr[])
{
  unsigned long i,ir=n,j,k,l=1;
  int jstack=0,*istack;
  double a,temp;
  
  istack=ivector(1,NSTACK);
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	for (i=j-1;i>=1;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	}
	arr[i+1]=a;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
	if (arr[l+1] > arr[ir]) {
	  SWAP(arr[l+1],arr[ir])
	    }
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	  }
      if (arr[l+1] > arr[l]) {
	SWAP(arr[l+1],arr[l])
	  }
      i=l+1;
      j=ir;
      a=arr[l];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j]);
      }
      arr[l]=arr[j];
      arr[j]=a;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort2(unsigned long n, double arr[], double brr[])
{
  unsigned long i,ir=n,j,k,l=1;
  int *istack,jstack=0;
  double a,b,temp;
  
  istack=ivector(1,NSTACK);
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=1;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (!jstack) {
	free_ivector(istack,1,NSTACK);
	return;
      }
      ir=istack[jstack];
      l=istack[jstack-1];
      jstack -= 2;
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
	SWAP(brr[k],brr[l+1])
	if (arr[l+1] > arr[ir]) {
	  SWAP(arr[l+1],arr[ir])
	    SWAP(brr[l+1],brr[ir])
	    }
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	  SWAP(brr[l],brr[ir])
	  }
      if (arr[l+1] > arr[l]) {
	SWAP(arr[l+1],arr[l])
	  SWAP(brr[l+1],brr[l])
	  }
      i=l+1;
      j=ir;
      a=arr[l];
      b=brr[l];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
	  SWAP(brr[i],brr[j])
	  }
      arr[l]=arr[j];
      arr[j]=a;
      brr[l]=brr[j];
      brr[j]=b;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in sort2.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}
#undef M
#undef NSTACK
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
#define EPS1 0.001
#define EPS2 1.0e-8

double probks(double alam)
{
  int j;
  double a2,fac=2.0,sum=0.0,term,termbf=0.0;
  
  a2 = -2.0*alam*alam;
  for (j=1;j<=100;j++) {
    term=fac*exp(a2*j*j);
    sum += term;
    if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
    fac = -fac;
    termbf=fabs(term);
  }
  return 1.0;
}
#undef EPS1
#undef EPS2
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
void ksone(double data[], unsigned long n, double (*func)(double), double *d,
	   double *prob)
{
  double probks(double alam);
  void sort(unsigned long n, double arr[]);
  unsigned long j;
  double dt,en,ff,fn,fo=0.0;
  
  sort(n,data);
  en=n;
  *d=0.0;
  for (j=1;j<=n;j++) {
    fn=j/en;
    ff=(*func)(data[j]);
    dt=FMAX(fabs(fo-ff),fabs(fn-ff));
    if (dt > *d) *d=dt;
    fo=fn;
  }
  en=sqrt(en);
  *prob=probks((en+0.12+0.11/en)*(*d));
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
void kstwo(double data1[], unsigned long n1, double data2[], unsigned long n2,
	   double *d, double *prob)
{
  double probks(double alam);
  void sort(unsigned long n, double arr[]);
  unsigned long j1=1,j2=1;
  double d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;
  
  sort(n1,data1);
  sort(n2,data2);
  en1=n1;
  en2=n2;
  *d=0.0;
  while (j1 <= n1 && j2 <= n2) {
    if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=j1++/en1;
    if (d2 <= d1) fn2=j2++/en2;
    if ((dt=fabs(fn2-fn1)) > *d) *d=dt;
  }
  en=sqrt(en1*en2/(en1+en2));
  *prob=probks((en+0.12+0.11/en)*(*d));
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;
  
  dif=fabs(x-xa[1]);
  c=vector(1,n);
  d=vector(1,n);
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free_vector(d,1,n);
  free_vector(c,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double), double a, double b, int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

double qromb(double (*func)(double), double a, double b)
{
  void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
  double trapzd(double (*func)(double), double a, double b, int n);
  void nrerror(char error_text[]);
  double ss,dss;
  double s[JMAXP+1],h[JMAXP+1];
  int j;
  
  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=trapzd(func,a,b,j);
    if (j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) < EPS*fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }
  nrerror("Too many steps in routine qromb");
  return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/
#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
  double gammln(double xx);
  void nrerror(char error_text[]);
  int n;
  double sum,del,ap;
  
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) nrerror("x less than 0 in routine gser");
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    nrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/
double gammq(double a, double x)
{
  void gcf(double *gammcf, double a, double x, double *gln);
  void gser(double *gamser, double a, double x, double *gln);
  void nrerror(char error_text[]);
  double gamser,gammcf,gln;
  
  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */

/*============================================================================*/
#define EPS 1.0e-10

void chsone(double bins[], double ebins[], int nbins, int knstrn, double *df,
	    double *chsq, double *prob)
{
  double gammq(double a, double x);
  void nrerror(char error_text[]);
  int j;
  double temp;
  
  *df=nbins-knstrn;
  *chsq=0.0;
  for (j=1;j<=nbins;j++) {
    if (ebins[j] <= 0.0) nrerror("Bad expected number in chsone");
    temp=bins[j]-ebins[j];
    *chsq += temp*temp/ebins[j];
  }
  //printf("%e  %e\n", *df, *chsq);
  *prob=gammq(0.5*(*df),0.5*(*chsq));
}
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
void chstwo(double bins1[], double bins2[], int nbins, int knstrn, double *df,
	    double *chsq, double *prob)
{
  double gammq(double a, double x);
  int j;
  double temp, r, s;
  
  *df=nbins-knstrn;
  *chsq=0.0;
  r = s = 0.0;
  for (j=1;j<=nbins;j++) {
    r += bins1[j];
    s += bins2[j];
  }
  for (j=1;j<=nbins;j++)
    if (bins1[j] == 0.0 && bins2[j] == 0.0)
      --(*df);
    else {
      temp=sqrt(s/r)*bins1[j] - sqrt(r/s)*bins2[j];
      *chsq += temp*temp/(bins1[j]+bins2[j]);
    }
  /* falls nur ein einziges bin belegt ist und knstrn = 1 (gleiche Anzahl) */
  if (*df == 0.0) {
    nrwarning("degree of freedom == 0 -> set to 1!");
    *df = 1;
  }
  *prob=gammq(0.5*(*df),0.5*(*chsq));
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                 MODIFIED FUNCTIONS FOR SHMMS                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/*============================================================================*/
/* modified KS-Test: ascending data and weights expected (w[max] == 1)! */ 
void ksone_weighted_shmm(double data[], double weights[], unsigned long n,
			 double (*func)(smodel*, int, double), smodel *smo,
			 int state, double *d, double *prob)
{
  double probks(double alam);
  unsigned long j;
  double dt,en,ff,fn,fo=0.0;
  
  en=n;
  *d=0.0;
  for (j=1;j<=n;j++) {
    /*** new */
    fn=weights[j];//fn=j/en;
    /***/
    ff=(*func)(smo, state, data[j]);
    dt=m_max(fabs(fo-ff),fabs(fn-ff));
    if (dt > *d) *d=dt;
    fo=fn;
  }
  en=sqrt(en);
  *prob=probks((en+0.12+0.11/en)*(*d));
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
#define FUNC(smo,s,c,x) ((*func)(smo,s,c,x))

double trapzd_cdfstate_shmm(double (*func)(smodel*, int, double, double), 
			    smodel* smo, int state, double c, double a, 
			    double b, int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=0.5*(b-a)*(FUNC(smo,state,c,a)+FUNC(smo,state,c,b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(smo,state,c,x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */


/*============================================================================*/
#define EPS 1.0e-6
#define JMAX 30
#define JMAXP (JMAX+1)
#define K 5

double qromb_cdfstate_shmm(double (*func)(smodel*,int,double,double), 
			   smodel *smo, int state, double c, double a, double b)
{
  void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
  double trapzd_cdfstate_shmm(double (*func)(smodel*, int, double, double),
			      smodel *smo, int state, double c,double a, 
			      double b, int n);
  void nrerror(char error_text[]);
  double ss,dss;
  double s[JMAXP+1],h[JMAXP+1];
  int j;
  
  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=trapzd_cdfstate_shmm(func,smo,state,c,a,b,j);
    if (j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= EPS*fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }
  printf("a = %e, b = %e, fabs(dss) = %e, EPS*fabs(ss) = %e, fabs(ss) = %e\n",
	 a,b, fabs(dss), EPS*fabs(ss), fabs(ss));
  nrerror("Too many steps in routine qromb");
  return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL|'`2. */
