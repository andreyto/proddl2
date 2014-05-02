//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef __USER_CST_H__
#define __USER_CST_H__


  //AT

/***********************************************************************
* double cost_function
*	This is the users cost function to optimize
*	(find the minimum).
*	cost_flag is set to TRUE if the parameter set
*	does not violates any constraints
*       parameter_lower_bound and parameter_upper_bound may be
*       adaptively changed during the search.
***********************************************************************/

//pointer to actual cost function - must be set from the outside:
extern double
(*pcost_function) (double *x,
	       double *parameter_lower_bound,
	       double *parameter_upper_bound,
	       double *cost_tangents,
	       double *cost_curvature,
	       ALLOC_INT * parameter_dimension,
	       int *parameter_int_real,
	       int *cost_flag,
	       int *exit_code,
	       USER_DEFINES * USER_OPTIONS);

#endif //__USER_CST_H__
