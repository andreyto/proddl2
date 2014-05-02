//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_GEOM_QUATFIT_CALL_H__
#define AT_GEOM_QUATFIT_CALL_H__

// Drivers to quatfit.hpp for some typical input data

#include <blitz/array.h>
#include <blitz/tinyvec.h>

#include "PRODDL/Common/iterator.hpp"
#include "PRODDL/Common/function.hpp"

#include "PRODDL/Geom/quatfit.hpp"

#include "PRODDL/Geom/geom_except.hpp"

namespace PRODDL { namespace Geom { namespace Align {

template<typename T_num, int n_dim>
typename QuatFitTraits<T_num>::Matr3 
quatFitArrays(blitz::Array<blitz::TinyVector<T_num,3>,n_dim> points1,
	      blitz::Array<blitz::TinyVector<T_num,3>,n_dim> points2) {

  if( points1.size() != points2.size() ) throw quatFit_error("quatFitArrays(points1,points2): points1.size() != points2.size()");
  PRODDL::gen_const<T_num> equal_weight(T_num(1));
  return quatFit<T_num>(points1.begin(),
			points1.end(),
			points2.begin(),
		 PRODDL::make_ipipe_iterator(equal_weight));
}

template<typename T_num, int n_dim,int n_dim_weights>
typename QuatFitTraits<T_num>::Matr3 
quatFitArrays(blitz::Array<blitz::TinyVector<T_num,3>,n_dim> points1,
	      blitz::Array<blitz::TinyVector<T_num,3>,n_dim> points2,
	      blitz::Array<T_num,n_dim_weights> weights) {

  if( points1.size() != points2.size() ) throw quatFit_error("quatFitArrays(points1,points2,weights): points1.size() != points2.size()");
  if( points1.size() != weights.size() ) throw quatFit_error("quatFitArrays(points1,points2,weights): points1.size() != weights.size()");
  return quatFit<T_num>(points1.begin(),
			points1.end(),
			points2.begin(),
			weights.begin());
}

} } } // namespace PRODDL { namespace Geom::Align

#endif // AT_GEOM_QUATFIT_CALL_H__

