//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_MATH_FFTW_TRAITS_H__
#define PRODDL_MATH_FFTW_TRAITS_H__

// FftwPlan class to switch at compile time between float and double versions
// of FFTW API

#include <fftw3.h>

#include <boost/shared_ptr.hpp>

#include "PRODDL/Common/logger.hpp"

namespace PRODDL {

  template<typename T_num>
  class FftwPlan;

  template<>
  class FftwPlan<float> {

  public:

    typedef float T_real;
    typedef fftwf_complex T_complex;

    typedef fftwf_plan Plan;

    struct PlanHolder {

      Plan plan;

      PlanHolder(const Plan _plan) : plan(_plan)
      {

	dbg::trace t1(DBG_HERE);

      }

      ~PlanHolder() {

	dbg::trace t1(DBG_HERE);

	fftwf_destroy_plan(plan);

      }

      Plan&
      get() {

	dbg::trace t1(DBG_HERE);

	return plan;

      }

      const
      Plan&
      get() const {

	dbg::trace t1(DBG_HERE);

	return plan;

      }

    }; // PlanHolder

    typedef boost::shared_ptr<PlanHolder> SharedPtrPlan;


    void
    dft_r2c(int rank, const int* n,
	    T_real* in,T_complex* out,
	    unsigned flags) {

      dbg::trace t1(DBG_HERE);

      ptrPlan = SharedPtrPlan(new PlanHolder(fftwf_plan_dft_r2c(rank,n,in,out,flags)));

    }
	
    void
    dft_c2r(int rank, const int* n,
	    T_complex* in,T_real* out,
	    unsigned flags) {

      dbg::trace t1(DBG_HERE);

      ptrPlan = SharedPtrPlan(new PlanHolder(fftwf_plan_dft_c2r(rank,n,in,out,flags)));

    }    

    void
    execute() {

      dbg::trace t1(DBG_HERE);

      fftwf_execute(ptrPlan->get());

    }

    Plan getPlan() const {

      dbg::trace t1(DBG_HERE);

      return ptrPlan->get();

    }

  protected:

    SharedPtrPlan ptrPlan;

  };

  template<>
  class FftwPlan<double> {

  public:

    typedef double T_real;
    typedef fftw_complex T_complex;

    typedef fftw_plan Plan;

    struct PlanHolder {

      Plan plan;

      PlanHolder(const Plan _plan) : plan(_plan)
      {}

      ~PlanHolder() {

	fftw_destroy_plan(plan);

      }

      Plan&
      get() {

	return plan;

      }

      const
      Plan&
      get() const {

	return plan;

      }

    }; // PlanHolder

    typedef boost::shared_ptr<PlanHolder> SharedPtrPlan;


    void
    dft_r2c(int rank, const int* n,
	    T_real* in,T_complex* out,
	    unsigned flags) {

      ptrPlan = SharedPtrPlan(new PlanHolder(fftw_plan_dft_r2c(rank,n,in,out,flags)));

    }
	
    void
    dft_c2r(int rank, const int* n,
	    T_complex* in,T_real* out,
	    unsigned flags) {

      ptrPlan = SharedPtrPlan(new PlanHolder(fftw_plan_dft_c2r(rank,n,in,out,flags)));

    }    

    void
    execute() {

      fftw_execute(ptrPlan->get());

    }

    Plan getPlan() const {

      return ptrPlan->get();

    }

  protected:

    SharedPtrPlan ptrPlan;

  };

}; // namespace PRODDL

#endif // PRODDL_MATH_FFTW_TRAITS_H__
