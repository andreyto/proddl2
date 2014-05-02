//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_OPTIM_MAC_H__
#define AT_OPTIM_MAC_H__

// interface class to David MacKay's optimizer

#include "PRODDL/Optim/macopt.hpp"

#include "PRODDL/Optim/optim_except.hpp"

#include "PRODDL/Optim/objective_func.hpp"

// for 'is_same<Type1,Type2>'
#include <boost/type_traits.hpp>

// for 'noncopyable'
#include <boost/utility.hpp>

// for 'copy'
#include <algorithm>

namespace PRODDL { namespace Optim {


  class OptimConjGradMac : protected macopt::Macopt, boost::noncopyable
  {

  public:
    typedef double T_num;
    typedef ObjectiveFuncD1<double> F;

  protected:

    // object that calculates gradient of function to be optimized
    // value of function is never needed in this optimizer

    F *f;

    //number of parameters

    int np;

    // redefine virtual function inhereted from Macopt:
    // 'func' finds objective function value;
    // 'dfunc' finds gradient.
    // Note that 'func' is not used within optimization process -
    // it is only called from 'Macopt::maccheckgrad()' that one
    // can optionally use to test that f is in agreement with
    // it's gradient.
  
    //callback functions used by optimization routine:

//     virtual double func(double* _p)
//     {
//       // Macopt uses subscripts with offset 1,
//       // we assume that F uses subscripts with offset 0

//       e = f->func(_p+1,np);
//       return e;
//     }

    virtual void dfunc(double* _p, double* _g)
    {

      // Macopt uses subscripts with offset 1,
      // we assume that F uses subscripts with offset 0

      //calculate objective function derivatives:

      f->grad(_p+1,_g+1,np);
    }
   
  public:

    OptimConjGradMac(F& f_,int _npar,double tol_=0.000001,int iter_max=1000,int _verb_flag=0) :
      Macopt(_npar,_verb_flag,tol_,iter_max),
      f(&f_), np(_npar)
    {
    }

    template<typename T_num> void optimize(T_num *p, int _np)
    {

      if( np != _np) throw optim_param_error("OptimConjGradmac::optimize(): Number of parameters does not match constructor");
      
      const bool is_double_arg = boost::is_same<T_num,double>::value;

      double *pbuf = 0;
      if( ! is_double_arg ) {
	pbuf = new double[np];
	std::copy(p,p+np,pbuf);
      }
      else {
	// cast is safe since we will get here only when T_num is double
	// cast is needed only to let this compile when T_num is not double
	pbuf = reinterpret_cast<double*>(p);
      }

      // optimize
      macoptII(pbuf-1,np);

      if( ! is_double_arg ) {
	std::copy(pbuf,pbuf+np,p);
	delete[] pbuf;
      }

    }

    inline int getSize() { return np; }  

    template<typename T_num>
    inline
    void setTolerance(T_num tol) {
      a_tol = double(tol);
    }

    inline
    void setIterMax(int itmax) {
      a_itmax = itmax;
    }

  };




} } // namespace PRODDL::Optim

#endif // AT_OPTIM_MAC_H__

  
