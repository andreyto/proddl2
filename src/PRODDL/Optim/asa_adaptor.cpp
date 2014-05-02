//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include <asa_user.h>

double f_test_asa(double *x,
	       double *parameter_lower_bound,
	       double *parameter_upper_bound,
	       double *cost_tangents,
	       double *cost_curvature,
	       ALLOC_INT * parameter_dimension,
	       int *parameter_int_real,
	       int *cost_flag,
	       int *exit_code,
	       USER_DEFINES * USER_OPTIONS)
{

  double q_n, d_i, s_i, t_i, z_i, c_r;
  int k_i;
  LONG_INT i, j;
  static LONG_INT funevals = 0;

  s_i = 0.2;
  t_i = 0.05;
  c_r = 0.15;

  q_n = 0.0;
  for (i = 0; i < *parameter_dimension; ++i)
    {
      j = i % 4;
      switch (j)
	{
	case 0:
	  d_i = 1.0;
	  break;
	case 1:
	  d_i = 1000.0;
	  break;
	case 2:
	  d_i = 10.0;
	  break;
	default:
	  d_i = 100.0;
	}
      if (x[i] > 0.0)
	{
	  k_i = (int) (x[i] / s_i + 0.5);
	}
      else if (x[i] < 0.0)
	{
	  k_i = (int) (x[i] / s_i - 0.5);
	}
      else
	{
	  k_i = 0;
	}

      if (fabs (k_i * s_i - x[i]) < t_i)
	{
	  if (k_i < 0)
	    {
	      z_i = k_i * s_i + t_i;
	    }
	  else if (k_i > 0)
	    {
	      z_i = k_i * s_i - t_i;
	    }
	  else
	    {
	      z_i = 0.0;
	    }
	  q_n += c_r * d_i * z_i * z_i;
	}
      else
	{
	  q_n += d_i * x[i] * x[i];
	}
    }
  funevals = funevals + 1;

  *cost_flag = TRUE;
  return (q_n);
}

#if 0

double f_test_asa_big(double *x,
	       double *parameter_lower_bound,
	       double *parameter_upper_bound,
	       double *cost_tangents,
	       double *cost_curvature,
	       ALLOC_INT * parameter_dimension,
	       int *parameter_int_real,
	       int *cost_flag,
	       int *exit_code,
	       USER_DEFINES * USER_OPTIONS)
  {
  /* Objective function from
   * %A A. Corana
   * %A M. Marchesi
   * %A C. Martini
   * %A S. Ridella
   * %T Minimizing multimodal functions of continuous variables
   *    with the "simulated annealing" algorithm
   * %J ACM Trans. Mathl. Software
   * %V 13
   * %N 3
   * %P 262-279
   * %D 1987
   *
   * This function, when used with ASA_TEST_POINT set to TRUE, contains
   * 1.0E20 local minima.  When *parameter_dimension is equal to 4, visiting
   * each minimum for a millisecond would take about the present age of the
   * universe to visit all these minima. */

  /* defines for the test problem, which assume *parameter_dimension
     is a multiple of 4.  If this is set to a large number, you
     likely should set Curvature_0 to TRUE. */
  double q_n, d_i, s_i, t_i, z_i, c_r;
  int k_i;
#if ASA_TEST_POINT
  int k_flag;
#endif
  LONG_INT i, j;
#if SELF_OPTIMIZE
#else
  static LONG_INT funevals = 0;
#endif
#if ASA_TEMPLATE_SAVE
  static int read_test = 0;
  FILE *ptr_read_test;
#endif

#if MY_TEMPLATE			/* MY_TEMPLATE_diminishing_ranges */
  /* insert code to automate changing ranges of parameters */
#endif
#if ASA_TEMPLATE		/* example of diminishing ranges */
  if (USER_OPTIONS->Locate_Cost == 12 && *(USER_OPTIONS->Best_Cost) < 1.0)
    {
      fprintf (ptr_out, "best_cost = %g\n", *(USER_OPTIONS->Best_Cost));
      for (i = 0; i < *parameter_dimension; ++i)
	{
	  parameter_lower_bound[i] = USER_OPTIONS->Best_Parameters[i]
	    - 0.5 * fabs (parameter_lower_bound[i]
			  - USER_OPTIONS->Best_Parameters[i]);
	  parameter_upper_bound[i] = USER_OPTIONS->Best_Parameters[i]
	    + 0.5 * fabs (parameter_upper_bound[i]
			  - USER_OPTIONS->Best_Parameters[i]);
	  parameter_lower_bound[i] = MIN (parameter_lower_bound[i],
				   USER_OPTIONS->Best_Parameters[i] - 0.01);
	  parameter_upper_bound[i] = MAX (parameter_upper_bound[i],
				   USER_OPTIONS->Best_Parameters[i] + 0.01);
	}
    }
#endif /* ASA_TEMPLATE */

  /* a_i = parameter_upper_bound[i] */
  s_i = 0.2;
  t_i = 0.05;
  c_r = 0.15;

#if ASA_TEST_POINT
  k_flag = 0;
  for (i = 0; i < *parameter_dimension; ++i)
    {
      if (x[i] > 0.0)
	{
	  k_i = (int) (x[i] / s_i + 0.5);
	}
      else if (x[i] < 0.0)
	{
	  k_i = (int) (x[i] / s_i - 0.5);
	}
      else
	{
	  k_i = 0;
	}
      if (k_i == 0)
	++k_flag;
    }
#endif /* ASA_TEST_POINT */

  q_n = 0.0;
  for (i = 0; i < *parameter_dimension; ++i)
    {
      j = i % 4;
      switch (j)
	{
	case 0:
	  d_i = 1.0;
	  break;
	case 1:
	  d_i = 1000.0;
	  break;
	case 2:
	  d_i = 10.0;
	  break;
	default:
	  d_i = 100.0;
	}
      if (x[i] > 0.0)
	{
	  k_i = (int) (x[i] / s_i + 0.5);
	}
      else if (x[i] < 0.0)
	{
	  k_i = (int) (x[i] / s_i - 0.5);
	}
      else
	{
	  k_i = 0;
	}

#if ASA_TEST_POINT
      if (fabs (k_i * s_i - x[i]) < t_i && k_flag != *parameter_dimension)
#else
      if (fabs (k_i * s_i - x[i]) < t_i)
#endif
	{
	  if (k_i < 0)
	    {
	      z_i = k_i * s_i + t_i;
	    }
	  else if (k_i > 0)
	    {
	      z_i = k_i * s_i - t_i;
	    }
	  else
	    {
	      z_i = 0.0;
	    }
	  q_n += c_r * d_i * z_i * z_i;
	}
      else
	{
	  q_n += d_i * x[i] * x[i];
	}
    }
  funevals = funevals + 1;

#if ASA_TEMPLATE_SAVE
  /* cause a crash */
  if ((ptr_read_test = fopen ("asa_save", "r")) == NULL)
    {
      read_test = 1;
      fclose (ptr_read_test);
    }
  else
    {
      fclose (ptr_read_test);
    }
  /* will need a few hundred if testing ASA_PARALLEL to get an asa_save */
  if (funevals == 50 && read_test == 1)
    {
      fprintf (ptr_out, "\n\n*** intended crash to test ASA_SAVE *** \n\n");
      fflush (ptr_out);
      printf ("\n\n*** intended crash to test ASA_SAVE *** \n\n");
      exit (2);
    }
#endif

  *cost_flag = TRUE;

#if SELF_OPTIMIZE
#else
#if TIME_CALC
  /* print the time every PRINT_FREQUENCY evaluations */
  if ((PRINT_FREQUENCY > 0) && ((funevals % PRINT_FREQUENCY) == 0))
    {
      fprintf (ptr_out, "funevals = %ld  ", funevals);
      print_time ("", ptr_out);
    }
#endif
#endif

#if ASA_TEMPLATE_SAMPLE
  USER_OPTIONS->Cost_Acceptance_Flag = TRUE;
  if (USER_OPTIONS->User_Acceptance_Flag == FALSE && *cost_flag == TRUE)
    USER_OPTIONS->Acceptance_Test (q_n, *parameter_dimension, USER_OPTIONS);
#endif /* ASA_TEMPLATE_SAMPLE */

  return (q_n);

#if ASA_TEMPLATE_SAMPLE

  int n;
  double cost;

  if (*cost_flag == FALSE)
    {
      for (n = 0; n < *parameter_dimension; ++n)
	cost_tangents[n] = 2.0 * x[n];
    }

  cost = 0.0;
  for (n = 0; n < *parameter_dimension; ++n)
    {
      cost += (x[n] * x[n]);
    }

  *cost_flag = TRUE;

  USER_OPTIONS->Cost_Acceptance_Flag = TRUE;
  if (USER_OPTIONS->User_Acceptance_Flag == FALSE && *cost_flag == TRUE)
    USER_OPTIONS->Acceptance_Test (cost, *parameter_dimension, USER_OPTIONS);

  return (cost);
#endif /* ASA_TEMPLATE_SAMPLE */
#if MY_TEMPLATE			/* MY_TEMPLATE_cost */
  /* Use the parameter values x[] and define your cost_function.
     The {} brackets around this function are already in place. */
#endif /* MY_TEMPLATE cost */
}

#endif //0
