//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_BZ_UNIFORM_ROT_H__
# define AT_BZ_UNIFORM_ROT_H__

/* 
   A.Tovchigrechko (2001)
   Blitz interface to "PRODDL/Geom/uniform_rot.hpp" functions to produce uniformly
   distributed rotations
*/

#include "PRODDL/Geom/uniform_rot.hpp"

#include <blitz/tinymat.h>
#include <blitz/tinyvec.h>
#include <blitz/array.h>


namespace PRODDL { namespace Geom {  namespace Transforms {

  // Return (3 x 3) rotation matrix built based on position 'x' in 3d rotational space.
  // See "PRODDL/Geom/uniform_rot.hpp" for detailed description.

  template<typename T_num> blitz::TinyMatrix<T_num,3,3> uniformRotation(const blitz::TinyVector<T_num,3>& x) {
    blitz::TinyMatrix<T_num,3,3> mat;
    uniformRotation(x.data(),mat.data());
    return mat;
  }


  // Return (N x 3 x 3) blitz::Array where each 3x3 slice is rotation matrix such that
  // matrices fill the 3d rotational space uniformly with step 'step'. Step must be in [0,1].

  template<typename T_num> blitz::Array<T_num,3> uniformRotationGridBz(T_num step) {
    int nM = uniformRotationGrid(step,(T_num*)0,0);
    blitz::Array<T_num,3> M(nM,3,3);
    uniformRotationGrid(step,M.data(),nM);
    return M;
  }
  
} } } // namespace PRODDL::Geom::Transforms

#endif // AT_BZ_UNIFORM_ROT_H__
