//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_OPTIM_OBJECTIVE_FUNC_BZ_H__
# define AT_OPTIM_OBJECTIVE_FUNC_BZ_H__

/*
  A.Tovchigrechko (2001)
  ObjectiveFunction adaptors that accept from optimizer flat T_num arrays and call
  the actual ('core') object that calculates the function with blitz::Array arguments
*/

#include "PRODDL/Optim/objective_func.hpp"

#include "PRODDL/Common/common_types.hpp"

namespace PRODDL { namespace Optim {

  class bz_c_cast_error : public std::exception {
  protected:
    std::string m_msg;
  public:
    bz_c_cast_error(const std::string& msg) throw(): 
      m_msg(msg) {
    }
    virtual const char* what() const throw() {
      return m_msg.c_str();
    }
    virtual ~bz_c_cast_error() throw (){}
  };

  // accepts T_num arrays;
  // casts them to blitz::Array<T_el,1> and calls
  // Ret (CoreType::*pmf_grad)(const blitz::Array<T_el,1>& points, blitz::Array<T_el,1>& gradients)
  // Requirements: T_el must be a wrapper object around a contiguous block of T_num elements, so that
  // sizeof(T_num)*(number of T_num elements in T_el) == sizeof(T_el). Also, accessing T_el as array
  // of T_num must be valid. Example of T_el is, of course, TinyVector<T_num,n_el>.


  template<typename T_num, class T_el, class CoreType,class Ret> class ObjectiveFuncContigBzArr : 
    public ObjectiveFuncD1<T_num> {

  public:

    typedef blitz::Array<T_el,1> T_array;
    typedef Ret (CoreType::*PmfGrad)(const T_array&, T_array&);

  protected:
    CoreType* m_data;
    PmfGrad m_pmf_grad;

  public:
    ObjectiveFuncContigBzArr(CoreType& data, PmfGrad pmf_grad) :
      m_data(&data),
      m_pmf_grad(pmf_grad) {
    }

    virtual void grad(const T_num *p, T_num *g, int np) {
      T_array bz_p = bz_cast(p,np);
      T_array bz_g = bz_cast(g,np);
      (m_data->*m_pmf_grad)(bz_p,bz_g);
    }

    static T_array bz_cast(const T_num *p, int np) {
      int n_el = np/(sizeof(T_el)/sizeof(T_num));
      return  T_array(const_cast<T_el*>(reinterpret_cast<const T_el*>(p)),blitz::shape(n_el),blitz::neverDeleteData);
    }

    static T_num* num_cast(const T_array& p, int& np) {
      if( ! p.isStorageContiguous() ) 
	throw bz_c_cast_error("ObjectiveFuncContigBzArr::num_cast(): Contiguous blitz Array is required to cast to C array");
      np = p.size()*(sizeof(T_el)/sizeof(T_num));
      return  const_cast<T_num*>(reinterpret_cast<const T_num*>(p.dataFirst()));
    }

  };


} } // namespace PRODDL::Optim

#endif // AT_OPTIM_OBJECTIVE_FUNC_BZ_H__
