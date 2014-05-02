//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_UNIFORM_POINTS_SPHERE3_H__
# define AT_UNIFORM_POINTS_SPHERE3_H__

/* 
   A.Tovchigrechko (2001)
   Generate points uniformly distributed on N-D sphere.
   Method 'uniformPointOnSphereOptim()' works this way:
     1. generate initial random set of N points inside a cube with side 1
     2. define interaction potential that consists of two parts:
       - spring to the center with equilibrium distance 1
       - repulsion to all other points
     2. minimize the interaction potential
*/


#include "PRODDL/Geom/geom_except.hpp"
#include <string>

#include <cmath>

#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/numinquire.h>
#include <random/uniform.h>

#include "PRODDL/Common/nd_index_iter.hpp"

#include "PRODDL/Common/bz_cast.hpp"

#include "PRODDL/Geom/uniform_points_sphere2.hpp"

#include "PRODDL/Optim/objective_func.hpp"

#include "PRODDL/Optim/optim_mac.hpp"

namespace PRODDL { namespace Geom { namespace Points {


  template<int n_dim,typename T_num>  class UniformPointsOnSphereOptim : public boost::noncopyable,
									 public PRODDL::Optim::ObjectiveFuncD1<double> {

  public:

    const T_num eps;
    const T_num eps_p2;
    const T_num dist_eps_p2;

    typedef blitz::TinyVector<T_num,n_dim> Point;

    typedef blitz::Array<T_num,2> APoints; // array N x n_dim

    typedef blitz::Array<Point,1> Points;

    typedef UniformPointsOnSphereOptim<n_dim,T_num> Self;

  protected:



    // points to distribute on a sphere
    Points m_points;

  public:


    UniformPointsOnSphereOptim(int n_points):
      eps(blitz::epsilon(T_num())),
      eps_p2(eps*eps),
      dist_eps_p2(eps_p2*4) {
    
      if( n_points <= 0 ) throw sphere_points_param_error<int>("uniformPointsOnSphereOptim(): n_points must be > 0","n_points",n_points);

      m_points.reference(Points(n_points));

    }

    // Generate initial set of points such that
    // no two points are closer than sqrt(dist_eps_p2).
    // It is needed to avoid having too large gradients at the beginning of
    // minimization process, or (worse) points with same coordinates,
    // in which later case they will stay together during minimization.

    void generateInitialPoints() {

      ranlib::UniformOpenClosed<T_num> rndgen;

      int n_points = m_points.rows();

      for(int i = 0; i < n_points; i++) {
	Point p_new;
	// We will try to generate new point up to 'max_precond_loops' times
	// till all existing points are no closer than sqrt(dist_eps_p2) to the new one.
	// If limit on number of loops is exceeded, exception 'sphere_points_loops_error' will
	// be thrown. We estimate that for any reasonable number of points (<1e5), probability
	// to get two random points close together is small ( total number of possible points
	// is ~ 1/sqrt(dist_eps_p2), which for 'float' is ~1e7; we assume that theoretical "granularity" of
	// random number generator is less than that of T_num - that appears to be true
	// for generator from Blitz, which has a period of 2^19937 - 1)
	const int max_precond_loops = 100;
	int loop_count = 0;
	for(bool too_close = true; too_close; loop_count++) {
	  too_close = false;
	  p_new = randomPointOnSphereCube<n_dim,T_num>(rndgen);
	  for(int j = i - 1; j >= 0; j--) {
	    Point& p_j = m_points(j);
	    Point d = p_new - p_j;
	    T_num dist_p2 = blitz::dot(d,d);
	    if(dist_p2 < dist_eps_p2) {
	      too_close = true;
	      break;
	    } // if(dist_p2
	  } // for(j
	  if(loop_count >= max_precond_loops)
	    throw PRODDL::Geom::sphere_points_loops_error("uniformPointsOnSphereOptim(): too many loops w/o result","max_precond_loops",loop_count);
	} // for(too_close
	m_points(i) = p_new;
      } // for(i
    } // generateInitialPoints()


    // compute potential's gradients 'grads' at position 'points'

    void getGrads(const Points& points, Points& grads) {

      grads = 0;
    
      // compute the spring-to-the-center component of potential
      // f(x) = k*(norm(x)-1)^2; grad f(x) = 2k*(x - (x/norm(x)) = 2k*x*(1 - (1/norm(x)));
      // we assume that x never gets too close to 0

      int n_points = points.size();

      T_num k_spring_twice = 2.0;
    
      for(int i=0; i < n_points; i++) {
	const Point& p = points(i);
	T_num r = std::sqrt(dot(p,p));
	grads(i) += ( k_spring_twice * (1 - (1/r)) ) * p;
      }

      // compute the mutual-repulsion component of the potential

      T_num k_repul = blitz::pow2(T_num(n_points));
      T_num r_cut_p2 = 16.0;
      if(n_points > 16) r_cut_p2 /= T_num(n_points);

      for(int i=0; i < n_points; i++) {
	const Point& p_i = points(i);
	Point& grad_i = grads(i);
	for(int j = i+1; j < n_points; j++) {
	  Point dif = points(j) - p_i;
	  T_num r_dif_p3 = dot(dif,dif);
	  if(r_dif_p3 < r_cut_p2) {
	    r_dif_p3 *= std::sqrt(r_dif_p3);
	    dif /= (r_dif_p3 * k_repul);
	    grad_i += dif;
	    grads(j) -= dif;
	  }
	}
      }

    } // getGrads()

    void optimize() {
      const T_num tolerance = 1e-6; //?eps*100;
      APoints aap(blitz::viewWithUnfoldedComponent<double>(m_points,1,0));
      PRODDL::Optim::OptimConjGradMac optimizer(*this,aap.size(),tolerance,1000,2);
      optimizer.optimize(aap.dataFirst(),aap.size());
    }

    void initAndOptimize() {
      generateInitialPoints();
      optimize();
    }

    virtual void grad(const double *p, double *g, int np) {
      APoints ap = APoints(const_cast<double*>(p),blitz::shape(np/n_dim,n_dim),blitz::neverDeleteData);
      APoints ag = APoints(                    g, blitz::shape(np/n_dim,n_dim),blitz::neverDeleteData);
      Points a_points_g(blitz::viewWithFoldedComponent<Point>(ag));
      getGrads(blitz::viewWithFoldedComponent<Point>(ap),a_points_g);
    }

    Points& getPoints() {
      return m_points;
    }

    const Points& getPoints() const {
      return m_points;
    }

  }; // class UniformPointsOnSphereOptim

  template<int n_dim,typename T_num> void uniformPointsOnSphereOptim(blitz::Array<blitz::TinyVector<T_num,n_dim>,1>& points) {
    UniformPointsOnSphereOptim<n_dim,T_num> points_gen(points.size());
    points_gen.initAndOptimize();
    points = points_gen.getPoints();
  }

  // Run checks on the set of points 'points' that are supposed to be uniformly distributed on
  // a sphere. Results: 'r_min' - must be supplied as array of same length as 'points'. It will 
  // contain distance to the closest other point for each point upon function completion.

  template<int n_dim,typename T_num> void checkUniformPointsOnSphere(const blitz::Array<blitz::TinyVector<T_num,n_dim>,1>& points,
								     blitz::Array<T_num,1>& r_min ) {

    typedef blitz::TinyVector<T_num,n_dim> Point;

    typedef blitz::Array<Point,1> Points;

    typedef blitz::Array<T_num,1> ANum;

    int n_p = points.size();

    // distance to the closest other point for each point

    if(n_p != r_min.size())
      throw PRODDL::Geom::sphere_points_error("checkUniformPointsOnSphere(): 'points' must have same length as 'r_min'");

    r_min = 0;

    for(int i = 0; i < n_p; i++) {
      const Point& p_i = points(i);
      for(int j = 0; j < i; j++) {
	const Point& p_j = points(j);
	Point d = p_i - p_j;
	T_num r2_ij = blitz::dot(d,d);
	if( r_min(i) > r2_ij ) r_min(i) = r2_ij;
	if( r_min(j) > r2_ij ) r_min(j) = r2_ij;
      }
    }

    r_min = blitz::sqrt(r_min);

  }

} } } // namespace PRODDL { namespace Geom::Points

#endif // AT_UNIFORM_POINTS_SPHERE3_H__
