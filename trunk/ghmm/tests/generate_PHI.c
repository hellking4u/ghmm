#include <ghmm/randvar.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>

/* bisection algorithm to explore precision*/
double biggest_thats_different_from_1()
{
  double low,up,half;
  low=0;
  up=100;
  while (up-low>0.001)
    {
      half=(low+up)/2.0;
      if (randvar_get_PHI(half)<1.0)
	low=half;
      else
	up=half;
    }
  printf("%f: biggest int, for which is PHI different from 1.0\n",low);
  return 0;
}

int print_PHI_table()
{
  long i;
  double a;
  gsl_sf_result r;
  i=0;
  for(a=0.0;a>-20.0;a-=0.001)
    {
      /* with gsl special function */
      (void)gsl_sf_erf_e(a,&r);
      /* with randvar function */
      printf("%f,%e,%e\n",a,r.val/2.0+0.5,r.err/2.0);
      i++;
    }
  printf("Erzeugte %d Werte\n",i);
}

int main()
{
  (void)biggest_thats_different_from_1();
  return 0;
}
