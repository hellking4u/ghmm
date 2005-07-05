#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ghmm/smodel.h>
#include <check.h>


/*******************************************************************************

  Constants for all test cases.

*******************************************************************************/

#define kSMODEL_SAMPLE_SMO_FILE "sample.smo"


/****** smodel_suite *********************************************************/
START_TEST (smodel_read_free) /* Add smodel_read_free in  smodel_suite() below */
{
  /* Test whether we can read an smodel, iterate over it, and free it again */

  smodel **model_array;
  int model_counter, result;
  char *inFileName = kSMODEL_SAMPLE_SMO_FILE;
  model_array = smodel_read(inFileName, &model_counter);

  fail_unless(model_counter == 1, 
	      "Read %d number of models instead of one", model_counter);

  fail_unless(model_array != NULL,
	      "smodel_read returned null pointer");

  result = smodel_free(model_array);  
  
  fail_unless(result == NULL,
	      "smodel_free failed");  
}
END_TEST







Suite *smodel_suite(void)
{
  Suite *s = suite_create("smodel");
  TCase *tc_core = tcase_create("Core");

  suite_add_tcase (s, tc_core); 
  
  tcase_add_test(tc_core, smodel_read_free);

  return s;
}
/****** end of smodel_suite ************************************************/




int main(void)
{
  int nf;
  Suite *s = smodel_suite();
  SRunner *sr = srunner_create(s);

  ghmm_rng_init();  /* Important! initialise rng  XXX UNTESTED*/

  srunner_run_all(sr, CK_NORMAL);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (nf == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

