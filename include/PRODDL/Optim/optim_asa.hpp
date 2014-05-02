//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_ASAOPT_H__
#define PRODDL_ASAOPT_H__

#include "PRODDL/Optim/optim_asa_base.hpp"

#include "PRODDL/Common/common_types.hpp"

// Concrete implementation of the class interface to ASA
// parametrized by the floating point type

namespace PRODDL { namespace Optim {

template<T_num> class Asa : public AsaBase {

public:

  typedef typename PRODDL::common_types::num_vector_type<T_num>::Type tVect;

  typedef ObjectiveFunc<T_num> ObjFuncType;

public:

  virtual double cost_function (double *x,
				double *parameter_lower_bound,
				double *parameter_upper_bound,
				double *cost_tangents,
				double *cost_curvature,
				ALLOC_INT * parameter_dimension,
				int *parameter_int_real,
				int *cost_flag,
				int *exit_code,
				USER_DEFINES * USER_OPTIONS) = 0;
  

public:

  T_num optimize(tVect& param,tVect& lowerBound,tVect& upperBound) {

    int n_param = param.size();

    asa_param.resize(n_param);

    asa_lower_bound.resize(n_param);

    asa_upper_bound.resize(n_param);

    asa_param = param;

    asa_lower_bound = lower_bound;

    asa_upper_bound = upper_bound;

    asa_exit_code = -1;

    asa_f = 0;
    
    globalAsaBase = this;

    asa_main(&asa_f,asa_param.data(),&asa_exit_code,asa_lower_bound.data(),
	     asa_upper_bound.data(),asa_param.size());

    param = asa_param;

    return T_num(asa_f);

  }

};

}} // namespace PRODDL::Optim

#endif // PRODDL_ASAOPT_H__
