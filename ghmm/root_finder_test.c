#include "root_finder.h"
#include <stdio.h>

/* find simple square root */
double test_sqrt(double x1, double x2, double x3, double x4)
{
  return x1*x1-x2;
}

double test_cubic(double x1, double x2, double x3, double x4)
{
  return x1*x1*x1+x1*x1*x2+x1*x3+x4;
}

int main()
{
  double result=0;

  result=zbrent_AB(test_sqrt,
		   0,10,0.000001,
		   2,0,0);

  printf("root 1: %f ?= 1.414214 \n",result);

  /* (x-2)^2*(x+5) */
  result=zbrent_AB(test_cubic,
		   -10,0,0.000001,
		   1,-16,20);

  printf("root 2: %f ?= -5\n",result);

  /* (x-1)*(x-2)*(x-3) */
  result=zbrent_AB(test_cubic,
		   0,1.5,0.000001,
		   -6,11,-6);

  printf("root 3: %f ?= 1\n",result);

  result=zbrent_AB(test_cubic,
		   1.5,2.5,0.000001,
		   -6,11,-6);

  printf("root 4: %f ?= 2\n",result);

  result=zbrent_AB(test_cubic,
		   2.5,5,0.000001,
		   -6,11,-6);

  printf("root 3: %f ?= 3\n",result);

  return 0;
}

