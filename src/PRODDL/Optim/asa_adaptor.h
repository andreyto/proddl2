//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef __ASA_ADAPTOR_H__
#define __ASA_ADAPTOR_H__

#include "user.h"
#include "user_cst.h"
#include <blitz/array.h>
#include <limits.h>
#include "mdebug.h"

using namespace blitz;

//basic adaptor class for ASA optimization package.
//uses static members initialized on each call to asa_run() method
//in order to supply asa_main() with callback function pointer.
//All members are protected (included constructor) so that
//this class to be used only as a base for some other class
//that defines the size of parameter array, provides the
//interface and T argument. Note that one can write:
//class child;
//class child : public asa_adaptor<child> {
//double f(...)
//child : asa_adaptor<child>(this) {
//   asa_param.resize(...);
//}
//void run() {
//   asa_run();
//}
//etc...
//};

//guard against some of asa defines that our interface does not support:
#if ! (MY_TEMPLATE && (!ASA_TEST) && ( OPTIONS_FILE) && ASA_LIB && ASA_TEMPLATE_LIB)
#error UNSUPPORTED COMBINATION OF ASA DEFINES, CHECK asa_user.h
#endif

template<class T> class asa_adaptor
{
protected:
typedef asa_adaptor Self;
typedef T cost_obj;
typedef Array<double,1> dVect;
  //parameters to be optimized by asa:
  dVect asa_param;
  dVect asa_lower_bound;
  dVect asa_upper_bound;
  int asa_exit_code;
  double asa_f;
  //set to this at each call to run(), so that we can call
  //pass address of a static callback function stat_cost_function
  //to asa_main():
  static cost_obj *p_stat_adaptor;
  cost_obj *p_adaptor;
  //address of this function is assigned to global pcost_function,
  //which is in turn used by asa:
  static double
  stat_cost_function (double *x,
	       double *parameter_lower_bound,
	       double *parameter_upper_bound,
	       double *cost_tangents,
	       double *cost_curvature,
	       ALLOC_INT * parameter_dimension,
	       int *parameter_int_real,
	       int *cost_flag,
	       int *exit_code,
	       USER_DEFINES * USER_OPTIONS)
{
  return Self::p_stat_adaptor->f(x,
				 parameter_lower_bound,
				 parameter_upper_bound,
				 cost_tangents,
				 cost_curvature,
				 parameter_dimension,
				 parameter_int_real,
				 cost_flag,
				 exit_code,
				 USER_OPTIONS);
}


  asa_adaptor(cost_obj *_p_adaptor) :
    p_adaptor(_p_adaptor), asa_exit_code(INT_MIN)
    {}
  //resets static and global variables and calls asa_main,
  //saves the results of optimization.
  void run_asa() {
    p_stat_adaptor = Self::p_adaptor;
    pcost_function = &Self::stat_cost_function;
    asa_main(&asa_f,asa_param.data(),&asa_exit_code,asa_lower_bound.data(),
	     asa_upper_bound.data(),asa_param.size());
    //DEBUG:
    //printf("\nIN MAIN FUNCTION:\n");
    //printf("exit_code=%d\n",asa_exit_code);
    //printf("f=%f\n",asa_f);
    //for(int i = 0; i < asa_param.size(); i++)
    //printf("fparam(%d)=%f\n",i,asa_param(i));
  }

public:
  
int exit_code_asa() const { return asa_exit_code; }
double cost_asa() const { return asa_f; }

//run asa with the TEST cost function from the ASA distribution's user_cst.h  as a regression test
  void test_asa() {
    //declare prototype - defintion is in asa_adaptor.cpp
    extern double f_test_asa(double *x,
	       double *parameter_lower_bound,
	       double *parameter_upper_bound,
	       double *cost_tangents,
	       double *cost_curvature,
	       ALLOC_INT * parameter_dimension,
	       int *parameter_int_real,
	       int *cost_flag,
	       int *exit_code,
	       USER_DEFINES * USER_OPTIONS);
    pcost_function = &f_test_asa;
    const int n_param = 4;
    dVect test_asa_param(n_param);
    dVect test_asa_lower_bound(n_param);
    dVect test_asa_upper_bound(n_param);
    test_asa_param = 999.0, -1007.0, 1001.0, 903.0;
    test_asa_lower_bound = -10000;
    test_asa_upper_bound =  10000;
    double test_asa_f=0;
    int test_asa_exit_code=0;
    asa_main(&test_asa_f,test_asa_param.data(),&test_asa_exit_code,test_asa_lower_bound.data(),
	     test_asa_upper_bound.data(),test_asa_param.size());
    OUTVAR(test_asa_f); OUTVAR(test_asa_param); OUTVAR(test_asa_exit_code); OUTENDL;
  }    

};

template<class T>  T *asa_adaptor<T>::p_stat_adaptor = 0;

#endif //__ASA_ADAPTOR_H__
