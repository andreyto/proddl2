//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Optim/macopt.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
   
   This version with some modifications by Mark Gibbs. (MNG)
*/

namespace macopt {

Macopt::Macopt(int n, 
	       int _verbose,
	       double _tolerance,
	       int _itmax,
	       int _rich) :
  a_n(n),
  a_tol(_tolerance),
  a_itmax(_itmax),
  a_verbose(_verbose),
  a_rich(_rich)
{
  a_g = new double[n+1] ;
  a_h = new double[n+1] ;
  a_xi = new double[n+1] ;
  a_pt = new double[n+1] ; /* scratch vector for sole use of macprod */
  a_gx = new double[n+1] ; /* scratch gradients             */
  a_gy = new double[n+1] ; /* used by maclinmin and macprod */
  /* if verbose = 1 then there is one report for each
     line minimization.
     if verbose = 2 then there is an additional 
     report for
     each step of the line minimization.
     if verbose = 3 then extra debugging 
     routines kick in.
     */

  a_end_if_small_step = 1 ; /* Change this to 0/1 if you prefer */
  a_stepmax = 0.5; // 0.5

  a_grad_tol_tiny = 1e-16 ; /* Probably not worth fiddling with */
  a_step_tol_tiny = 0.0 ;   /* Probably not worth fiddling with */
  a_linmin_maxits = 20 ;    /* Probably not worth fiddling with */
  a_lastx = 0.01 ;          /* only has a transient effect      */
  a_lastx_default = 0.01 ;  /* -- defines typical distance in parameter
				space at which the line minimum is expected;
				both these should be set. the default is 
				consulted if something goes badly wrong and 
				a reset is demanded. */

/* don't fiddle with the following, unless you really mean it */
  a_linmin_g1 = 2.0 ; 
  a_linmin_g2 = 1.25 ; 
  a_linmin_g3 = 0.5 ; 
  a_restart = 0 ; 
}

Macopt::~Macopt() {
  delete[] a_g;
  delete[] a_h;
  delete[] a_xi;
  delete[] a_pt;
  delete[] a_gx;
  delete[] a_gy;
}


void Macopt::macoptII
  (double *p,            /* starting vector                                */
   int    n             /* number of dimensions                           */
   )                   
{
  int j ;
  double gg , gam , dgg ;
  double *g , *h , *xi ;
  int end_if_small_grad = 1 - a_end_if_small_step ;
  double step , tmpd ;

  /* A total of 7 double * 1..n are used by this optimizer. 
     p           is provided when the optimizer is called 
     pt          is used by the line minimizer as the temporary vector. 
                    this could be cut out with minor rewriting, using p alone
     g, h and xi are used by the cg method as in NR - could one of these
                     be cut out?
     the line minimizer uses an extra gx and gy to evaluate two gradients. 
     */

  g = a_g ; 
  h = a_h ;
  xi = a_xi ;
  
  dfunc( p , xi );

  macopt_restart ( 1 ) ; 

  for ( a_its = 1 ; a_its <= a_itmax ; a_its ++ ) {

    for ( gg = 0.0 , j = 1 ; j <= n ; j ++ ) 
      gg += g[j]*g[j];          /* find the magnitude of the old gradient */
    a_gtyp = sqrt ( gg / (double)(n) ) ; 

    if ( a_verbose > 0 ) 
      printf ( "mac_it %d of %d : gg = %6.3g tol = %6.3g: ", a_its , a_itmax , gg , a_tol ) ;

    if ( ( end_if_small_grad && ( gg <= a_tol ) ) 
	|| ( gg <= a_grad_tol_tiny ) ) {
      //      macopt_free ( a ) ;
      if ( a_verbose > 0 ) printf ("\n");
      return;
    }

    step = maclinminII ( p  ) ; 

    //if (step > a_stepmax) {
    //  if ( a_verbose > 1 ) printf (" (step %9.5g truncated to stepmax)",step);
    //  step = a_stepmax;
    //}

    if ( a_restart == 0 ) {
      if ( a_verbose > 1 ) printf (" (step %9.5g)",step);
      if ( a_verbose > 0 ) printf ("\n");
      if ( ( a_end_if_small_step  && ( step <= a_tol ) ) 
	  || ( step <= a_step_tol_tiny ) ) {
	//	macopt_free ( a ) ;
	return;
      }
    }

    /* if we are feeling rich, evaluate the gradient at the new
       `minimum'. alternatively, linmin has already estimated this
       gradient by linear combination of the last two evaluations and
       left it in xi */
    if ( a_rich || a_restart ) { 
      dfunc( p , xi  ) ; 
    }
    if ( a_restart ) {
      if ( a_verbose > 0 ) fprintf(stderr,"Restarting macopt\n" ) ; 
      macopt_restart ( 0 ) ;
/* this is not quite right
   should distinguish whether there was an overrun indicating that the 
   value of lastx needs to be bigger / smaller; 
   in which case resetting lastx to default value may be a bad idea, 
   giving an endless loop of resets 
*/
    } else {
      dgg=0.0;
      for ( j = 1 ; j <= n ; j ++ ) {
	dgg += ( xi[j] + g[j] ) * xi[j] ;
      }
      gam = dgg / gg ;
      for ( tmpd = 0.0 , j = 1 ; j <= n ; j ++ ) {
	g[j] = -xi[j];                /* g stores (-) the most recent gradient */
	xi[j] = h[j] = g[j] + gam * h[j] ;
	/* h stores xi, the current line direction */
	/* check that the inner product of gradient and line search is < 0 */
	tmpd -= xi[j] * g[j] ; 
      }
      if ( tmpd > 0.0  || a_verbose > 2 ) {
	if ( a_verbose > 0 ) fprintf(stderr,"new line search has inner prod %9.4g\n", tmpd ) ; 
      }
      if ( tmpd > 0.0 ) { 
	if ( a_rich == 0 ) {
	  if ( a_verbose > 0 ) fprintf (stderr, "Setting rich to 1; " ) ; 
	  a_rich = 1 ; 
	}
	a_restart = 2 ; /* signifies that g[j] = -xi[j] is already done */
	if ( a_verbose > 0 ) fprintf(stderr,"Restarting macopt (2)\n" ) ; 
	macopt_restart ( 0 ) ;
      }
    }
  }
  if ( a_verbose > 0 ) fprintf(stderr,"Reached iteration limit in macopt; continuing.\n"); 
  //  macopt_free ( a ) ;	
  return;
} /* NB this leaves the best value of p in the p vector, but
     the function has not been evaluated there if rich=0     */

/* maclinmin.
   Method: 
       evaluate gradient at a sequence of points and calculate the inner 
       product with the line search direction. Continue until a 
       bracketing is achieved ( i.e a change in sign ). */
double Macopt::maclinminII 
(
 double *p 
 )
{
  int n = a_n ; 

  double x , y ;
  double s , t , m ;
  int    its = 1 , i ;
  double step , tmpd ; 
  double  *gx = a_gx , *gy = a_gy ;

  /* at x=0, the gradient (uphill) satisfies s < 0 */
  if ( a_verbose > 2 ) { /* check this is true: (no need to do this really
			    as it is already checked at the end of the main
			    loop of macopt) */

    tmpd = macprodII ( p , gy , 0.0 );

    if ( a_verbose > 0 ) {
      fprintf (stderr, "inner product at 0 = %9.4g\n" ,tmpd ) ;
    }
    if ( tmpd > 0.0 ) { 
      a_restart = 1 ; 
      return 0.0 ; 
    }
  }

  x = a_lastx / a_gtyp ;

  /* MNG alteration :

  Here we check to make sure that the step in x is not going
  to change the value of the hyperparameters by anything more
  than 1.0 on the first go. This is to avoid initial massive
  steps that are making things screw up
  */
  
  double max = 0.0;
  double x_new = x;
  for(i=1; i <= n; ++i){
    double temp1 = fabs(a_xi[i]*x);
    if( temp1 > max && temp1 > 1.0){
      max = temp1;
      x_new = fabs(1.0 / a_xi[i]);
    }
  }
  if( x != x_new ){
    if(a_verbose > 1)
      fprintf(stderr,"Note : maclinmin step truncated\n");
    x = x_new;
  }

  s = macprodII ( p , gx , x ) ; 
  
  if ( s < 0 )  {  /* we need to go further */
    do {
      y = x * a_linmin_g1 ;
      t = macprodII ( p , gy , y  ) ; 
      if ( a_verbose > 1 ) 
	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
      if ( t >= 0.0 ) break ;
      x = y ; s = t ; a_gunused = gx ; gx = gy ; gy = a_gunused ; 
      its++ ;
      /* replaces: for ( i = 1 ; i <= n ; i ++ ) gx[i] = gy[i] ; */
    }
    while ( its <= a_linmin_maxits ) ;
  } else if ( s > 0 ) { /* need to step back inside interval */
    do {
      y = x * a_linmin_g3 ;
      t = macprodII ( p , gy , y ) ; 
      if ( a_verbose > 1 ) 
	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
      if ( t <= 0.0 ) break ;
      x = y ; s = t ; a_gunused = gx ; gx = gy ; gy = a_gunused ; 
      its ++ ;
    } while ( its <= a_linmin_maxits ) ;
  } else { /* hole in one s = 0.0 */
    t = 1.0 ; y = x;
  }

  if ( its > a_linmin_maxits )  {
    if ( a_verbose > 0 ) fprintf (stderr, "Warning! maclinmin overran" );
    /* this can happen where the function goes \_ and doesn't buck up
       again; it also happens if the initial `gradient' does not satisfy
       gradient.`gradient' > 0, so that there is no minimum in the supposed
       downhill direction.  I don't know if this actually happens... If it
       does then I guess a_rich should be 1.

       If the overrun is because too big a step was taken then
       the interpolation should be made between zero and the most 
       recent measurement. 

       If the overrun is because too small a step was taken then 
       the best place to go is the most distant point. 
       I will assume that this doesn't happen for the moment.

       Also need to check up what happens to t and s in the case of overrun.
       And gx and gy. 

       Maybe sort this out when writing a macopt that makes use of the gradient 
       at zero? 
    */
    tmpd = macprodII ( p , gy , 0.0 );
    if ( a_verbose > 0 ) {
      fprintf (stderr, "- inner product at 0 = %9.4g\n" ,
	       tmpd ) ; 
    }
    if ( tmpd > 0 && a_rich == 0 ) {
      if ( a_verbose > 0 ) fprintf (stderr, "setting rich to 1\n" ) ;       
      a_rich = 1 ; 
    }
    if ( tmpd > 0 ) a_restart = 1 ; 
  }


  /*  Linear interpolate between the last two. 
      This assumes that x and y do bracket. */
  if ( s < 0.0 ) s = - s ;
  if ( t < 0.0 ) t = - t ;
  m = ( s + t ) ;
  s /= m ; t /= m ;
  
  m =  s * y + t * x ; 
  /* evaluate the step length, not that it necessarily means anything */
  for ( step = 0.0 , i = 1 ; i <= n ; i ++ ) {
    tmpd = m * a_xi[i] ;
    p[i] += tmpd ; /* this is the point where the parameter vector steps */
    step += fabs ( tmpd ) ; 
    a_xi[i] = s * gy[i] + t * gx[i] ;
    /* send back the estimated gradient in xi (NB not like linmin) */
  }
  a_lastx = m * a_linmin_g2 *  a_gtyp ;
  
  return ( step / (double) ( n ) ) ; 
}

double Macopt::macprodII 
( 
 double *p , double *gy , double y  

 ) {
  double *pt = a_pt ; 
  double *xi = a_xi ; 
  /* finds pt = p + y xi and gets gy there, 
     returning gy . xi */
  int n = a_n ; 

  int i;
  double s = 0.0 ;

  for ( i = 1 ; i <= n ; i ++ ) 
    pt[i] = p[i] + y * xi[i] ;
  
  dfunc( pt , gy ) ;

  for ( i = 1 ; i <= n ; i ++ ) 
    s += gy[i] * xi[i] ;

  return s ;
}

void Macopt::macopt_restart ( int start ) 
/* if start == 1 then this is the start of a fresh macopt, not a restart */
{
  int j , n=a_n ; 
  double *g, *h, *xi ;
  g = a_g ;  h = a_h ;  xi = a_xi ; 

  if ( start == 0 ) a_lastx = a_lastx_default ; 
  /* it is assumed that dfunc( p , xi  ) ;  has happened */
  for ( j = 1 ; j <= n ; j ++ ) {
    if ( a_restart != 2 ) g[j] = -xi[j] ;
    xi[j] = h[j] = g[j] ;
  }
  a_restart = 0 ; 
}

} // namespace macopt

/*
<!-- hhmts start -->
Last modified: Tue Dec 17 17:49:44 1996
<!-- hhmts end -->
*/

