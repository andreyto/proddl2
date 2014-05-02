#ifndef AT_NR_ARR_H__
#define AT_NR_ARR_H__

// Adaptor class to expose 2D array to NRC (numerical recipes in C) routines.
// NRC requires "pointer to array of pointers" for 2D array.
// This class takes pointer to flat memory area and array bounds, and constructs array of pointers to 
// rows.

#include <boost/utility.hpp>

namespace PRODDL { namespace nr {

  template<typename T_num> class MatrixAdaptor : boost::noncopyable
  {

  public:
    typedef T_num T_number;

  protected:

    T_num **a; //array of row pointers:a=0-crl, a[i]=0-ccl

    int crl,crh,ccl,cch; //low and high indices, bounds inclusive

  public:

    inline int rl() {  return crl;  }
    inline int rh() {  return crh;  }
    inline int cl() {  return ccl;  }
    inline int ch() {  return cch;  }
    inline int rows() {  return (crh-crl+1);  }
    inline int cols() {  return (cch-ccl+1);  }
 

    // 'adataFirst' - points to first data element, [crl_,crh_] - row index range (closed),
    // [ccl_,cch_] - column index range (closed)

    MatrixAdaptor(const T_num* adataFirst, int crl_,int crh_,int ccl_,int cch_):
      crl(crl_), crh(crh_), ccl(ccl_), cch(cch_)
    {
      int crows=(crh-crl+1);
      int ccols=(cch-ccl+1);
      a = new T_num*[crows];
      T_num* p = const_cast<T_num*>(adataFirst) - ccl;
      for(int i=0; i<crows; i++, p+=ccols)
	{
	  a[i]=p;
	}
      a-=crl;
    }

    ~MatrixAdaptor() {
      delete[] (a+crl);
    }

    // return array of pointers to rows, so that individual element can be accessed as
    // rowPointers()[i][j], where i is in [crl_,crh_], j is in [ccl_,cch_].

    T_num** rowPointers() {
      return a;
    }
	
    const T_num** rowPointers() const {
      return a;
    }
  };

} } // namespace PRODDL::nr

#endif // AT_NR_ARR_H__
