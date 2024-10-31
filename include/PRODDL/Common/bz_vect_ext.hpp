//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_BZ_VECT_EXT_H__
#define AT_BZ_VECT_EXT_H__

// Some extensions for handling vectors in space in Blitz:
// e.g. norm of vector

#include <limits> // Blitz still include limits.h but uses std::
#include <blitz/array.h>
#include <blitz/tinyvec2.h>

#include <cmath>

namespace blitz_ext {

  // Power 2 of Euclidian norm

  template<class T_expr> inline typename T_expr::T_numtype
  dotSelf(const T_expr& v)
  {
    //return blitz::dot(v,v);
    return blitz::sum(v*v);
  }

  // Euclidian norm

  template<class T_expr> inline typename T_expr::T_numtype
  normSelf(const T_expr& v)
  {
    return std::sqrt(blitz::dot(v,v));
  }


} // namespace blitz_ext

#endif // AT_BZ_VECT_EXT_H__
