//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_UNIFORM_POINTS_SPHERE2_H__
# define AT_UNIFORM_POINTS_SPHERE2_H__

/* 
   A.Tovchigrechko (2001)
   Generate points uniformly distributed on 3D sphere.
   Method 'uniformPointOnSphereCube()' works this way:
     1. generate points uniformly inside a cube
     2. eliminate all points outside the unit sphere
     2. project remaining points on a unit sphere by rescaling their norm to 1
   Method 'randomPointOnSphereCube()' creates one random point on a spehere
   by generating random point taken from a uniform distribution
   inside a cube untill the point falls inside the sphere,
   then rescaling it to length 1
*/


#include "PRODDL/Geom/geom_except.hpp"
#include <string>

#include <math.h>


#include "PRODDL/Common/common_types.hpp"

#include <random/uniform.h>


#include "PRODDL/Common/nd_index_iter.hpp"

namespace PRODDL { namespace Geom { namespace Points {

  template<int n_dim,typename T_num> struct UniformPointsOnSphereCube {

    typedef blitz::TinyVector<int,n_dim> IndGrid;
    typedef blitz::TinyVector<T_num,n_dim> Coord;
    
    const int n_side;
    const T_num f_side, R, R_p2;
    const T_num f_step;

    UniformPointsOnSphereCube(int _n_side): 
      n_side(_n_side),
      f_side(1.0), R(f_side), R_p2(R*R), f_step(f_side/n_side) {
    }    

    int getNMaxPoints() const {
      int n = n_side*2 + 1;
      return n*n*n;
    }

    template<typename OutIter> OutIter operator()(OutIter res) const {
      if( n_side <= 0 ) throw sphere_points_error("uniformPointsOnSphereCube(): n_side must be > 0");
      IndGrid ind_grid, lbound(-n_side), ubound(n_side);
	  ubound += 1; // index_mover uses open range
      typedef PRODDL::index_mover<n_dim> ind_mover;
      for(ind_mover::before_first(ind_grid,lbound); ind_mover::next(ind_grid,lbound,ubound); ) {
	Coord point_grid = ind_grid * f_step;
	T_num r_p2 = blitz::dot(point_grid,point_grid);
	if(r_p2 <= R_p2 && r_p2 != 0) {
	  point_grid /= sqrt(r_p2);
	  *res = point_grid;
	  ++res;
	}
      }
      return res;      
    }

  }; // struct 


  // Return random n_dim point within cube of [-1,1] or (-1,1), with or w/o hole at (0), depending 
  // on the random number generator supplied.
  // Argument: random number generator of uniform numbers within (0,1) or (0,1] etc.

  template<int n_dim,typename T_num,class RndGen> blitz::TinyVector<T_num,n_dim> 
  randomPoint(RndGen& rndgen) {
    blitz::TinyVector<T_num,n_dim> point;
    for(int i=0; i < n_dim; i++) point(i) = (rndgen.random()-0.5)*2.0;
  }

  // Return a random point within a sphere of a given radius (default 1.0)
  // Point coords are taken from a random number generator 'rndgen',
  // which must return values within (0,1), (0,1] etc.
  // Values provided by 'rndgen' are scaled to fall into cube (-1,1)
  // and then subject to acceptance criterion that checks if
  // the point is actually within a sphere of radius 1.

  template<int n_dim,
	   typename T_num,
	   class RndGen> 
  blitz::TinyVector<T_num,n_dim> 
  randomPointInSphereCube(RndGen& rndgen, T_num radius = 1.0) {
    const T_num eps_p2 = blitz::pow2(blitz::epsilon(T_num()));
    blitz::TinyVector<T_num,n_dim> point;
    while(true) {
      for(int i=0; i < n_dim; i++) point(i) = (rndgen.random() - 0.5)*2.0;
      T_num r_p2 = blitz::dot(point,point);
      if( r_p2 <= 1 ) {
	point *= radius;
	break;
      }
    }
    return point;
  }

  // Output result of calling 'randomPointInSphereCube'
  // n times into 'out' iterator range.
  // RndGen template parameter must have a default ctor.

  template<int n_dim, typename T_num, class OutputIterator>
  OutputIterator
  randomPointsInSphereCube(T_num radius, int n, OutputIterator out) {
    ranlib::UniformOpenClosed<T_num> rndgen;
    for( int i = 0; i < n; ++i, ++out ) {
      *out = randomPointInSphereCube<n_dim>(rndgen,radius);
    }
    return out;
  }


  // Return random n_dim point on a sphere of radius 1.
  // Argument: random number generator of uniform numbers within (0,1) or (0,1] etc.

  template<int n_dim,typename T_num,class RndGen> blitz::TinyVector<T_num,n_dim> 
  randomPointOnSphereCube(RndGen& rndgen) {
    const T_num eps_p2 = blitz::pow2(blitz::epsilon(T_num()));
    blitz::TinyVector<T_num,n_dim> point;
    while(true) {
      for(int i=0; i < n_dim; i++) point(i) = (rndgen.random() - 0.5)*2.0;
      T_num r_p2 = blitz::dot(point,point);
      if( r_p2 <= 1 && r_p2 >= eps_p2 ) {
	point /= sqrt(r_p2);
	break;
      }
    }
    return point;
  }
  
} } } // namespace PRODDL { namespace Geom::Points

#endif // AT_UNIFORM_POINTS_SPHERE2_H__
