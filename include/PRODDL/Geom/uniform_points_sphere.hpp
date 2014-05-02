//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_UNIFORM_POINTS_SPHERE_H__
# define AT_UNIFORM_POINTS_SPHERE_H__

/* 
   A.Tovchigrechko (2001)
   Generate points uniformly distributed on 3D sphere.
   Method 'uniformPointOnSphereSpiral()' is 
   an adaptation of code from  Joseph O'Rourke.
   His description of method is retained below.

*/


#include "PRODDL/Geom/geom_except.hpp"
#include <string>

#include "PRODDL/Common/c_array.hpp"

#include <cmath>

namespace PRODDL { namespace Geom { namespace Points {


/*
This program will generate a given number of spiral points uniformly 
distributed on the surface of a sphere.
The idea behind the algorithm is that one can cut the globe with 
N horizontal planes spaced 2/(N-1) units apart, forming N circles of 
latitude on the sphere, each latitude containing one spiral point.  To
obtain the kth spiral point, one proceeds upward from the (k-1)st point
(theta(k-1), phi(k-1)) along a great circle to the next latitude and 
travels counterclockwise along ti for a fixed distance to arrive at the 
kth point (theta(k), phi(k)).

Reference: E.B. Saff and A.B.J. Kuijlaars, Distributing Many Points on a Sphere,
The Mathematical Intelligencer, 19(1), Winter (1997);

Written by Joseph O'Rourke and Min Xu, June 1997.
Used in the textbook, "Computational Geometry in C."
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/


  template<typename T_num> void uniformPointsOnSphereSpiral(T_num* outPoints, int nPoints) {

    if( nPoints > 0 ) {

      // cast 'ares' to T_num[][3]
      typename PRODDL::c_array_traits<T_num,2,3>::A_type aPoints = PRODDL::c_array_cast<3>(outPoints);

      T_num phi1, phi, theta, h, x, y, z;

      // AT:the original code claims that it will generate 'n' points, but obviously generates
      // only n-1. this is the blunt work-around.
      int n = nPoints+1;

      phi1 = 0.0;


      T_num *point = aPoints[0];

      for (int k = 2; k <= n - 1; k ++ ) {
	/* Generate a random point on a sphere of radius 1. */
	h = -1 + 2 * ( k - 1 ) / ( T_num )( n - 1 );
	theta = std::acos ( h );

	if ( theta < 0 || theta > M_PI )
	  throw sphere_points_error("uniformPointsOnSphereSpiral(): Unexpected error - computed theta out of allowed domain");

	phi = phi1 + 3.6 / ( std::sqrt ( ( T_num )n * ( 1 - h * h ) ) ); 
	phi = std::fmod ( phi, 2 * M_PI );
	phi1 = phi;

	x = std::cos ( phi ) * std::sin ( theta );
	y = std::sin ( phi ) * std::sin ( theta );
	/* z = std::cos ( theta ); But z==h, so: */
	z = h;

	point = aPoints[k-2];
	point[0] = x; point[1] = y; point[2] = z;

      } // for

      // last point
      point = aPoints[nPoints - 1];
      point[0] = 0; point[1] = 0; point[2] = 1;

    } // if (nPoints>0)

  } // function



} } } // namespace PRODDL { namespace Geom::Points

#endif // AT_UNIFORM_POINTS_SPHERE_H__
