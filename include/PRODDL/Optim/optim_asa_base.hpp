//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_ASAOPT_BASE_H__
#define PRODDL_ASAOPT_BASE_H__

#include "PRODDL/Common/common_types.hpp"

// Class interface to ASA

namespace PRODDL { namespace Optim {

class AsaBase {

protected:

  extern AsaBase * globalAsaBase;

  typedef PRODDL::common_types::num_vector_type<double>::Type dVect;

  //parameters to be optimized by asa:

  dVect asa_param;
  dVect asa_lower_bound;
  dVect asa_upper_bound;

  int asa_exit_code;

  double asa_f;

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
  

protected:

  void optimize() {
    
    globalAsaBase = this;

    asa_main(&asa_f,asa_param.data(),&asa_exit_code,asa_lower_bound.data(),
	     asa_upper_bound.data(),asa_param.size());
  }

  double getValue() {
    return asa_f;
  }

  double getExitCode() {
    return asa_exit_code;
  }

};

}} // namespace PRODDL::Optim

#endif // PRODDL_ASAOPT_BASE_H__
