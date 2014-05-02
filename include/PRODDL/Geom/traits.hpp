//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_GEOM_TRAITS_H__
#define AT_GEOM_TRAITS_H__

#include "PRODDL/Common/common_types.hpp"

namespace PRODDL { namespace Geom { 

  template<typename T_num> struct SpaceTraits {
    typedef typename PRODDL::common_types::small_matrix_type<T_num,3>::Type Matrix3x3;
    typedef typename PRODDL::common_types::point_type<T_num,3>::Type Point3;
    typedef typename PRODDL::common_types::point_type<int,3>::Type IntPoint3;
    typedef typename PRODDL::common_types::point_type<Point3,2>::Type Point3Pair;
    typedef typename PRODDL::common_types::point_type<Point3,3>::Type Point3Triple;
    typedef typename PRODDL::common_types::num_vector_of_points_type<T_num,3>::Type VPoint3;
    typedef typename PRODDL::common_types::num_multiarray_type<T_num,2>::Type LargeMatrix;
  };

} } //namespace PRODDL { namespace Geom {

#endif // AT_GEOM_TRAITS_H__
