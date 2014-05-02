//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_COMMON_TYPES_H__
#define AT_COMMON_TYPES_H__

// Type generating classes to define some commonly used types, such as
// points in space. The idea is that code must include this file,
// and do something like:
// typedef typename at::common_types::point_type<T_num,n_dim>::Type Point;
// And use 'Point' as point type. Thus, it will be possible to replace
// the undelying type in this file with minimal changes to other code, for
// instance, go from blitz::TinyVector to Pooma.

#include <blitz/tinyvec2.h>
#include <blitz/tinymat2.h>
#include <blitz/array.h>

#include <vector>

namespace PRODDL { namespace common_types {

  template<typename T_num,int n_dim> struct point_type {

    typedef blitz::TinyVector<T_num,n_dim> Type;

  };

  template<class T_num,int n_dim> struct small_vector_type {

    typedef blitz::TinyVector<T_num,n_dim> Type;

  };

  template<typename T_num,int n_rows, int n_cols=n_rows> struct small_matrix_type {

    typedef blitz::TinyMatrix<T_num,n_rows,n_cols> Type;

  };

  template<typename T_num> struct num_vector_type {

    typedef blitz::Array<T_num,1> Type;

  };

  template<typename T_num> struct num_matrix_type {

    typedef blitz::Array<T_num,2> Type;

  };

  struct num_index_type {

    typedef blitz::Array<int,1> Type;

  };

  // dynamic (growing) vector type. example is std::vector

  template<typename T_num> struct dyn_vector_type {

    typedef std::vector<T_num> Type;

  };

  template<typename T_num,int n_dim> struct num_multiarray_type {

    typedef blitz::Array<T_num,n_dim> Type;

  };

  template<typename T_num, int n_dim> struct num_vector_of_points_type {

    typedef typename num_vector_type<typename point_type<T_num,n_dim>::Type>::Type Type;

  };

} } // namespace PRODDL::common_types

#endif // AT_COMMON_TYPES_H__
