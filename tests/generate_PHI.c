#include <ghmm/randvar.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>

int main()
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

