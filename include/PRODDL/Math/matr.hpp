//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_MATH_MATR_H__
#define PRODDL_MATH_MATR_H__

// Some matrix manipulation routines

#include "PRODDL/Geom/traits.hpp"

#include "PRODDL/exceptions.hpp"

#include "PRODDL/Blitz/bztinyvec.hpp"

namespace PRODDL { namespace Math {

  // inject 'transpose()' for blitz::TinyMatrix
  using namespace blitz_ext;

  template<typename T_num> bool inverse (const typename PRODDL::Geom::SpaceTraits<T_num>::Matrix3x3& matrIn,
					 typename PRODDL::Geom::SpaceTraits<T_num>::Matrix3x3& matrOut, 
					 T_num tolerance = 1.e-6)
  {

    matrOut(0,0) = matrIn(1,1)*matrIn(2,2) -
      matrIn(1,2)*matrIn(2,1);
    matrOut(0,1) = matrIn(0,2)*matrIn(2,1) -
      matrIn(0,1)*matrIn(2,2);
    matrOut(0,2) = matrIn(0,1)*matrIn(1,2) -
      matrIn(0,2)*matrIn(1,1);
    matrOut(1,0) = matrIn(1,2)*matrIn(2,0) -
      matrIn(1,0)*matrIn(2,2);
    matrOut(1,1) = matrIn(0,0)*matrIn(2,2) -
      matrIn(0,2)*matrIn(2,0);
    matrOut(1,2) = matrIn(0,2)*matrIn(1,0) -
      matrIn(0,0)*matrIn(1,2);
    matrOut(2,0) = matrIn(1,0)*matrIn(2,1) -
      matrIn(1,1)*matrIn(2,0);
    matrOut(2,1) = matrIn(0,1)*matrIn(2,0) -
      matrIn(0,0)*matrIn(2,1);
    matrOut(2,2) = matrIn(0,0)*matrIn(1,1) -
      matrIn(0,1)*matrIn(1,0);

    T_num det =
      matrIn(0,0)*matrOut(0,0) +
      matrIn(0,1)*matrOut(1,0)+
      matrIn(0,2)*matrOut(2,0);

    if ( fabs(det) <= tolerance )
      return false;

    T_num invDet = 1.0/det;
    //TODO: a cludge until the proper operator support is available
    //for blitz::TinyMatr
    for(int i = 0; i < 3; i++) for(int j=0; j < 3; j++)
      matrOut(i,j) *= invDet;
    return true;
  }

  template<typename T_num> typename PRODDL::Geom::SpaceTraits<T_num>::Matrix3x3 
  inverse (const typename PRODDL::Geom::SpaceTraits<T_num>::Matrix3x3& matrIn,
	 T_num tolerance = 1.e-6)
  {
    typename PRODDL::Geom::SpaceTraits<T_num>::Matrix3x3 matrOut;
    if (! inverse(matrIn,matrOut,tolerance))
      throw tolerance_error("inverse(): Tolerance exceeded");
    return matrOut;
  }


}} // namespace PRODDL::Math

#endif //PRODDL_MATH_MATR_H__
