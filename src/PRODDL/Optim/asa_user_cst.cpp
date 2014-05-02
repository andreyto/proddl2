//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//


#include "asa_user.h"

double
(*pcost_function) (double *x,
	       double *parameter_lower_bound,
	       double *parameter_upper_bound,
	       double *cost_tangents,
	       double *cost_curvature,
	       ALLOC_INT * parameter_dimension,
	       int *parameter_int_real,
	       int *cost_flag,
	       int *exit_code,
	       USER_DEFINES * USER_OPTIONS) = 0;


//my implementation of predefined name cost_function()
double
cost_function (double *x,
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
  //just call the pointer, must be already set
  return pcost_function(x,
		 parameter_lower_bound,
	       parameter_upper_bound,
	       cost_tangents,
	       cost_curvature,
	       parameter_dimension,
	       parameter_int_real,
	       cost_flag,
	       exit_code,
	       USER_OPTIONS);
}
