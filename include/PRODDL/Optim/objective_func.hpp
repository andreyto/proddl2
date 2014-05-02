//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_OPTIM_OBJECTIVE_FUNC_H__
# define AT_OPTIM_OBJECTIVE_FUNC_H__

namespace PRODDL { namespace Optim {


  // Object that calculates objective function
  // This object is passed to optimization object

  template<typename _T_num> class ObjectiveFunc {
  public:
    typedef _T_num T_num;

    // return function value for T_num p[np] array of
    // parameters

    virtual T_num func(const T_num *p, int np) = 0;

    virtual ~ObjectiveFunc() {}
  };


  // Object that calculates first derivative of objective function
  // w/o calculating the function itself
  // This object is passed to optimization object (e.g. OptimConjGradMac)

  template<typename _T_num> class ObjectiveFuncD1 {
  public:
    typedef _T_num T_num;


    // accept T_num p[np] array of parameters and
    // store gradient into T_num g[np] array

    virtual void grad(const _T_num *p, _T_num *g, int np) = 0;

    virtual ~ObjectiveFuncD1() {}
  };


  template<typename _T_num, typename T_ret> class ObjectiveFuncD1FromPtr : public ObjectiveFuncD1<_T_num> {

  public:
    
    // type of function pointer
    typedef T_ret (*FPtr)(const _T_num*,_T_num*,int);

  protected:
    
    FPtr fptr;

  public:
    
    explicit ObjectiveFuncD1FromPtr(FPtr _fptr) :
      fptr(_fptr) {
    }

    // define virtual inhereted methods

    // call function for T_num p[np] array of
    // parameters and store gradient to T_num g[np] array

    virtual void grad(const _T_num *p, _T_num *g, int np) {
      fptr(p,g,np);
    }

  };

  // The following will not compile in gcc 2.96 when instantiated.
//   // Functions to make creation of ObjectiveFuncD1FromPtr easier

//   template<typename T_num, typename T_ret> ObjectiveFuncD1FromPtr<T_num,T_ret>
//   makeObjectiveFuncD1FromPtr( T_ret (*fptr)(const T_num*,T_num*,int) ) {
//     return ObjectiveFuncD1FromPtr<T_num,T_ret>(fptr);
//   }

} } // namespace PRODDL::Optim

#endif // AT_OPTIM_OBJECTIVE_FUNC_H__
