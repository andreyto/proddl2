//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_OPTIONS_H__
#define PRODDL_OPTIONS_H__

namespace PRODDL {

  // arrays of rank up to five are supported; default template parameters below are
  // just a trick to allow usage like 'c_array_traits<T,rank,n_2,...,n_rank>'.
  // Compiler is supposed to issue error if number of dimension size parameters specified is less than
  // rank-1, since each dimension size of C array is required to be positive.

  template<typename T,int rank,int n2, int n3=-1, int n4=-1, int n5=-1> struct c_array_traits;

  template<typename T,int n2> struct c_array_traits<T,2,n2>
  {
    typedef T (*A_type)[n2];    
  };

  template<typename T,int n2, int n3> struct c_array_traits<T,3,n2,n3>
  {
    typedef T (*A_type)[n2][n3];    
  };

  template<typename T,int n2, int n3, int n4> struct c_array_traits<T,4,n2,n3,n4>
  {
    typedef T (*A_type)[n2][n3][n4];    
  };

  template<typename T,int n2, int n3, int n4, int n5> struct c_array_traits<T,5,n2,n3,n4,n5>
  {
    typedef T (*A_type)[n2][n3][n4][n5];    
  };

  template<typename T,int rank,int n2, int n3, int n4, int n5> struct c_array_traits
  {
    typedef void A_type;
  };

  template<int n2,typename T> typename c_array_traits<T,2,n2>::A_type 
  c_array_cast(T* a) {
    return reinterpret_cast<typename c_array_traits<T,2,n2>::A_type>(a);
  }

  template<int n2,int n3,typename T> typename c_array_traits<T,3,n2,n3>::A_type 
  c_array_cast(T* a) {
    return reinterpret_cast<typename c_array_traits<T,3,n2,n3>::A_type>(a);
  }

  template<int n2,int n3,int n4,typename T> typename c_array_traits<T,4,n2,n3,n4>::A_type 
  c_array_cast(T* a) {
    return reinterpret_cast<typename c_array_traits<T,4,n2,n3,n4>::A_type>(a);
  }

  template<int n2,int n3,int n4,int n5,typename T> typename c_array_traits<T,5,n2,n3,n4,n5>::A_type 
  c_array_cast(T* a) {
    return reinterpret_cast<typename c_array_traits<T,5,n2,n3,n4,n5>::A_type>(a);
  }

} // namespace PRODDL

#endif // AT_C_ARRAY_H__
