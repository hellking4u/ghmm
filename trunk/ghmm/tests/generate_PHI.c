#include <math.h>
#include <gsl/gsl_sf.h>
#include <ghmm/rng.h>
#include <ghmm/randvar.h>
#include <stdio.h>
#include <stdlib.h>

/* helper function and structure for rng tests */

typedef struct {
  double start;
  double step;
  double end;
  size_t num;
  unsigned int* bins;
} bin_row;

/* I like object orientation! */

int clear_bin_row(bin_row* bins)
{
  int i;
  for (i=0; i<bins->num; i++) bins->bins[i]=0;
  return i;
}

int measure_bin_row(bin_row* bins, double start, double end)
{
      bins->start=start;
      bins->end=end;
      bins->step=(end-start)/bins->num;
}

bin_row* create_bin_row(double start, double end, size_t num)
{
  bin_row* bins=NULL;
  /* sanity check */
  if (start==end || num<=0) return NULL;
  bins=(bin_row*)malloc(sizeof(bin_row)+sizeof(int)*num);
  if (bins==NULL)
    return bins;
  else
    {
      bins->num=num;
      bins->bins=(unsigned int*)((void*)bins+sizeof(bin_row));
      (void)measure_bin_row(bins,start,end);
      (void)clear_bin_row(bins);
      return bins;
    }
}

void delete_bin_row(bin_row* bins)
{
  if (bins!=NULL) free(bins);
}


int count_in_bin_row(bin_row* bins, double value)
{
  if (value<bins->start || value>bins->end)
    return -1;
  else 
    {
      size_t bin_num;
      bin_num=floor((value-bins->start)/bins->step);
      /* paranoid for very small/big numbers */
      if (bin_num>=bins->num) return -1;
      bins->bins[bin_num]++;
      return 0;
    }
}

int lineprint_bin_row(bin_row* bins, FILE* file)
{
  int i;
  for (i=0; i<bins->num; i++)
    fprintf(file,"%f\t%d\n",i*bins->step+bins->start,bins->bins[i]);
  return i;
}

/* bisection algorithms to explore precision*/
int greatest_thats_different_from_1(void)
{
  double low,up,half;
  low=0;
  up=10;
  /* find start interval */
  while (randvar_get_PHI(up)<1.0)
    {
      low=up;
      up*=2;
    }
  /* shrink it */
  while (up-low>0.001)
    {
      half=(low+up)/2.0;
      if (randvar_get_PHI(half)<1.0)
	low=half;
      else
	up=half;
    }
  printf("%f: greatest int, for which PHI is different from 1.0\n",low);
  return 0;
}

int least_thats_different_from_0(void)
{
  double low,up,half;
  low=-10;
  up=0;
  /* find start interval */
  while (randvar_get_PHI(low)>0.0)
    {
      up=low;
      low*=2;
    }
  /* shrink it */
  while (up-low>0.001)
    {
      half=(low+up)/2.0;
      if (randvar_get_PHI(half)>0.0)
	up=half;
      else
	low=half;
    }
  printf("%f: least int, for which PHI is different from 0.0\n",low);
  return 0;
}


int print_PHI_table(void)
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

void test_bins(void)
{
  int i;
  bin_row* bins;
  bins=create_bin_row(-10,10,100);

  gsl_rng_init();
  fprintf(stdout,"plot '-'\n");
  for (i=0; i<100000; i++)
    count_in_bin_row(bins,randvar_std_normal(0));
  (void)lineprint_bin_row(bins,stdout);
  fprintf(stdout,"e\n");

  delete_bin_row(bins);
}


/* create a gnuplot file with pictures from randomnumbers*/
void gnuplot_show_distributions(void)
{
  /* where output goes */
  FILE* file;
  const int generate_count=10000;
  /* count how often generator was used*/
  int generate_counter;
  /* row of bins */
  bin_row* bins;

  /* now to stdout*/
  file=stdout;
  /* init bin_row */
  bins=create_bin_row(-10,10,100);
  /* initialise generator */
  gsl_rng_init();

  /* prepare output for gnuplot */
  fprintf(file,"set terminal png small color\n");

  /* generate a normal distribution picture */
  fprintf(file,"set output 'std_normal_dist.png'; plot '-' title 'standard normal distribution'\n");
  for (generate_counter=0; generate_counter<generate_count; generate_counter++)
    count_in_bin_row(bins,randvar_std_normal(0));
  (void)lineprint_bin_row(bins,file);
  fprintf(file,"e\n");

  /* generate  two truncated normal distribution pictures */
  measure_bin_row(bins,-10,10);
  fprintf(file,"set output 'trunc_normal_dist.png';\n\
plot '-' title '1st','-' title '2nd', '-' title '3rd'\n");

  clear_bin_row(bins);
  for (generate_counter=0; generate_counter<generate_count; generate_counter++)
    /*    count_in_bin_row(bins,randvar_normal_pos(0,1,0)); */
    count_in_bin_row(bins,gsl_ran_gaussian_tail(RNG,-4,1));
  (void)lineprint_bin_row(bins,file);
  fprintf(file,"e\n");

  clear_bin_row(bins);
  for (generate_counter=0; generate_counter<generate_count; generate_counter++)
    count_in_bin_row(bins,gsl_ran_gaussian_tail(RNG,0,1));
  (void)lineprint_bin_row(bins,file);
  fprintf(file,"e\n");

  clear_bin_row(bins);
  for (generate_counter=0; generate_counter<generate_count; generate_counter++)
    count_in_bin_row(bins,gsl_ran_gaussian_tail(RNG,4,1));
  (void)lineprint_bin_row(bins,file);
  fprintf(file,"e\n");

  /* terminate gnuplot */
  fprintf(file,"set output; quit\n");

  /* cleanup */
  delete_bin_row(bins);
  
}

int main()
{
#if 0
  (void)greatest_thats_different_from_1();
  (void)least_thats_different_from_0();
#endif

  gnuplot_show_distributions();
  return 0;
}
