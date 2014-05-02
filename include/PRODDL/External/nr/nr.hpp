#ifndef AT_NR_H__
# define AT_NR_H__

/* 
   Andrei Tovchigrechko (2001)
   Declarations of templatized versions of NR functions.
   Instantiations for different types (e.g. float, double) must be provided elsewhere.
*/


namespace PRODDL { namespace nr {

  template<typename T_num>
  void jacobi(T_num **a, int n, T_num d[], T_num **v, int *nrot);

  template<typename T_num>
  void eigsrt(T_num d[], T_num **v, int n, bool ascend=false);

  template<typename T_num>
  T_num betacf(T_num a, T_num b, T_num x);

  template<typename T_num>
  T_num betai(T_num a, T_num b, T_num x);

  template<typename T_num>
  T_num gammln(T_num xx);

  template<typename T_num>
  T_num zbrent(T_num (*func)(T_num), T_num x1, T_num x2, T_num tol);

} } // namespace PRODDL::nr

#endif // AT_NR_H__

