//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef GDIAM_SIMPLE_H_
# define GDIAM_SIMPLE_H_

/* Simplified, type-neutral interface to GDIAM code for bounding box*/

namespace PRODDL { namespace Geom { namespace Gdiam {

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
			   double eps=0.0);

// Function to compute a tight-fitting bounding box of the point-set
//
// Input Parameters:
// double* points - 3xN array:
// a11 a12 a13
// a21 a22 a23
// ...
// where each row is one 3D point
// 
// int n - number of points
// int grid_n - grid size (number elements of grid in each dimension, used to reduce
// pairwise distance computations similar to MD codes?)
// int sample_n - number of points to sample as a representative subset
//
// Output Parameters:
// double* left, double* right - ends of diagonal, double[3] each;
// double* directions - 3 normals of box vertixes, double[9] 

void simp_approx_mvbb_grid_sample(const double* points, int n, 
				  double* left, double* right,
				  double* directions,
				  int grid_n, int sample_n);

} } } // namespace PRODDL::Geom::Gdiam

#endif // GDIAM_SIMPLE_H_
