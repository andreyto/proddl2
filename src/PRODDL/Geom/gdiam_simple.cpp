//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Geom/gdiam_simple.hpp"

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <memory.h>
#include  <math.h>

#include  "PRODDL/External/Gdiam/gdiam.h"

namespace PRODDL { namespace Geom { namespace Gdiam {

/* Simplified interface to GDIAM code for bounding box*/

//
// Info on GDIAM types (from looking at "gdiam.h" and "gdiam.cpp" contents):
// 'gdiam_real' is 'double'
// 'points' is array of 'double' of 3xN length, if N is a number of points
//
// GPointPair   gdiam_approx_diam( gdiam_point  * start, int  size,
//                                gdiam_real  eps )
//
// gdiam_bbox   gdiam_approx_mvbb_grid_sample( gdiam_point  * start, int  size,
//                                            int  grid_size, int  sample_size )
//


// Function to compute the diameter of a set of points
//
// Input Parameters:
// double* points - 3xN array:
// a11 a12 a13
// a21 a22 a23
// ...
// where each row is one 3D point
// 
// int n - number of points
// double eps - epsilon (?)
//
// Output Parameters:
// double* left and double* right each must point to
// array of size 3. The function will store coords for
// the ends of "diameter" line

void simp_approx_diam_pair(const double* points, int n, double* left, double* right, 
		      double eps) {

    GPointPair   pair;  
    pair = gdiam_approx_diam_pair( const_cast<double*>(points), n, eps );
    for(int i=0; i < 3; i++) {
      left [i] = pair.p[i];
      right[i] = pair.q[i];
    }
}


inline double* copy3dpoint(const double *src, double *dst) {
  for(int i = 0; i < 3; i++) dst[i] = src[i];
  return dst;
}

inline double* copy3dpoint(double x1, double x2, double x3, double* dst) {
  dst[0] = x1;
  dst[1] = x2;
  dst[2] = x3;
  return dst;
}

// Function to compute a tight-fitting bounding box of the point-set
//
// Input Parameters:
// double* points - Nx3 array:
// a11 a12 a13
// a21 a22 a23
// ...
// where each row is one 3D point
// 
// int n - number of points
// int grid_n - grid size (number elements of grid in each dimension?)
// int sample_n - number of points to sample as a representative subset
//
// Output Parameters:
// double* left, double* right - ends of diagonal, double[3] each;
// double* directions - 3 normals of box vertixes, double[9] 

void simp_approx_mvbb_grid_sample(const double* points, int n, 
				  double* left, double* right,
				  double* directions,
				  int grid_n, int sample_n) {

    gdiam_bbox   bb;

    gdiam_point  * pnt_arr = gdiam_convert( const_cast<double*>(points), n );

    bb = gdiam_approx_mvbb_grid_sample( pnt_arr, n, grid_n, sample_n );

    bb.get_vertex(0,0,0,left);
    bb.get_vertex(1,1,1,right);

    //copy3dpoint(bb.low_1,bb.low_2,bb.low_3,left);
    //copy3dpoint(bb.high_1,bb.high_2,bb.high_3,right);
    for(int i = 0; i < 3; i++) {
      copy3dpoint(bb.get_dir(i),directions+i*3);
    }

    free(pnt_arr); // gdiam_convert uses malloc(...) to allocate pnt_arr
    //bb.dump();

}


} } } // namespace PRODDL { namespace Geom::Gdiam
