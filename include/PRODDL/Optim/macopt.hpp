//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_ORIG_MACOPT_H_
#define AT_ORIG_MACOPT_H_

/* 
   
   http://131.111.48.24/mackay/c/macopt.html       mackay@mrao.cam.ac.uk

   Please do not use macopt without understanding a little about how it works;
   there are some control parameters which the user MUST set!

   David MacKay's optimizer, based on conjugate gradient ideas, 
   but using bracketing of the zero of the inner product 

             (gradient).(line_search_direction)

   to do the line minimization. Only derivative calculations are required.
   The length of the first step in the line search (often set to "1.0"
   in other code) is adapted here so that, if 0.00001 is a better step size, 
   it soon cottons on to that and saves ~log(10000) bracketing operations.
   The result is that (with rich set to 0) the program can use 
   as few as 2 derivatives per line search. (If rich is set to 1, it does 
   an extra derivative calculation at the beginning of each line search 
   making a minimum of 3 per line search. Set rich=0 if you think 
   that the surface is locally quite quadratic.) If the program does average 
   2 derivatives per line search then it must be superior to most cg methods 
   including use of Rbackprop (which costs 2 derivatives straight off)

   A possible modification: where the function can be returned at same 
   time as the dfunction --- there is nothing clever to do with the 
   value, but it could be used as a sanity check and a convergence criterion. 

   See http://131.111.48.24/mackay/c/macopt.html for further discussion.

   NB: The value of "tol" is totally arbitrary and must be set by 
   you to a value that works well for your problem. 
   It depends completely on the typical value of the gradient / step size. 

   Tol specifies a magnitude of gradient at which a halt is called. 
   or a step size.

   This program MINIMIZES a function.

   **********************************************

   Modified and converted to C++ class by Steve Waterhouse 8th April
   1997.

   **********************************************
*/

namespace macopt
{

class Macopt 
{
public:

  Macopt(int _n, 
	 int _verbose = 0, 
	 double _tolerance = 0.001, 
	 int _itmax = 100,
	 int _rich = 1);
  virtual ~Macopt();
  void macoptII(double *p, int    dim ); 

protected:
  //virtual double func(double* _p) = 0;
  virtual void dfunc(double* _p, double* _xi)  = 0;
  int a_n;     /* dimension of parameter space */
  double a_tol ;    /* convergence declared when the gradient vector is smaller
		     in magnitude than this, or when the mean absolute 
		     step is less than this (see above) */

  int a_itmax ;             /* max iterations */

private:
  double maclinminII(double *p);
  double macprodII (double * , double * , double ) ;
  void macopt_restart ( int ) ;

private:

  int a_verbose ; 


  double a_grad_tol_tiny ; /* if gradient is less than this, we definitely 
			    stop, even if we are not using a gradient 
			    tolerance */
  double a_step_tol_tiny ; /* if step is less than this, we stop, even if 
			    we are not using a step tolerance */
  int a_end_if_small_step ; /* defines the role of tol -- alternative is
			     end_on_small_grad */
  int a_its ;               /* number of its */

  int a_rich ; /* whether to do the extra gradient evaluation at the beginning 
	      of each new line min */
  double a_stepmax ;        /* largest step permitted (not used in macopt) */

  int a_linmin_maxits ;     /* in maclinmin */
  double a_linmin_g1 ;      /* factors for growing and shrinking the interval */
  double a_linmin_g2 ;
  double a_linmin_g3 ;
  double a_lastx     ;      /* keeps track of typical step length */
  double a_lastx_default ;  /* if maclinmin is reset, lastx is set to this */

/* These should not be touched by the user. They are handy pointers for macopt
   to use 
*/
  double a_gtyp ; /* stores the rms gradient for linmin */
  double *a_pt , *a_gx , *a_gy , *a_gunused ;
  double *a_xi , *a_g , *a_h ; 

  int a_restart ;           /* whether to restart macopt - fresh cg directions */

/* a_lastx :--- 1.0 might make general sense, (cf N.R.)
				  but the best setting of all is to have 
				  a prior idea of the eigenvalues. If 
				  the objective function is equal to sum of N
				  terms then set this to 1/N, for example 
				  Err on the small side to be conservative. */
};	

} //namespace macopt

#endif // AT_ORIG_MACOPT_H_

