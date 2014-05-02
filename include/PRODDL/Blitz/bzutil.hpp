//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_BZUTIL_H__
#define AT_BZUTIL_H__

//some extensions to Blitz

#include <blitz/array.h>

#include <algorithm>

#include "PRODDL/Blitz/bzarray_iter.hpp"


namespace blitz_ext
{

  template<typename T_num,int n_rank>
  blitz::GeneralArrayStorage<n_rank> 
  getArrayStorage(const blitz::Array<T_num,n_rank>& arr)
  {


    typedef blitz::GeneralArrayStorage<n_rank> Storage;
    Storage storage;

    for(int i = 0; i < n_rank; i++) storage.ascendingFlag()(i) = arr.isRankStoredAscending(i);
    storage.ordering() = arr.ordering();
    storage.base() = arr.base();

    return storage;

  }


  template<typename T_num_out,typename T_num_inp,int n_rank>
  blitz::Array<T_num_out,n_rank>
  makeArrayWithSameStructure(const blitz::Array<T_num_inp,n_rank>& arr)
  {

    return blitz::Array<T_num_out,n_rank>(arr.shape(),getArrayStorage(arr));

  }

  //Assign Array<T_num_out,n_rank> = Array<T_num_inp,n_rank>
  //Blitz does not have the corresponding operator when T_num_out <> T_num_inp
  //This function relies on slow STD iterators. It requires that 
  //operator= (T_num_out,T_num_inp) is defined.

  template<typename T_num_inp,typename T_num_out,int n_rank>
  void
  copyArrayFirstToSecond(const blitz::Array<T_num_inp,n_rank>& arr_inp,
	      blitz::Array<T_num_out,n_rank>& arr_out)
  {

    std::copy(arr_inp.begin(),arr_inp.end(),arr_out.begin());

  }  

  // Return contiguous reference to array's data
  // If input array is already contiguous, the
  // returned array will just reference existing data,
  // otherwise the contiguous copy will be made.
  // Because it is not known in advance if the copy
  // or the true reference will be returned, the 
  // returned value should be used read only and not assumed
  // to reflect changes later made to the original array's data.
  // Requirements:
  // operator T_num_return = T_num_input must be defined (hint:
  // it is defined for TinyVector, not defined for Array

  template<typename T_num_return, typename T_num_input,int n_dim>
  struct ArrayConverter; 

  template<typename T_num_input,int n_dim> 
  struct ArrayConverter<T_num_input,T_num_input,n_dim> {

    static
    blitz::Array<T_num_input,n_dim>
    getContiguousRefOrCopy(const blitz::Array<T_num_input,n_dim>& x) {
      if( x.isStorageContiguous() )
	return x;
      else
	return x.copy();
    }

  };

  template<typename T_num_return, typename T_num_input,int n_dim>
  struct ArrayConverter {

    static
    blitz::Array<T_num_return,n_dim>
    getContiguousRefOrCopy(const blitz::Array<T_num_input,n_dim>& x) {
      blitz::Array<T_num_return,n_dim> ret(makeArrayWithSameStructure<T_num_return>(x));
      //NOTE: 'ret = x' would not work - blitz::Array has copy operator and copy ctor only
      //for the expressions with the same T_num. Therefore, we use slow std iterators here
      copyArrayFirstToSecond(x,ret);
      return ret;
    }

  };

  template<typename T_num_return, typename T_num_input,int n_dim>
  blitz::Array<T_num_return,n_dim>
  getContiguousRefOrCopy(const blitz::Array<T_num_input,n_dim>& x) {
    return ArrayConverter<T_num_return,T_num_input,n_dim>::getContiguousRefOrCopy(x);
  }

} //namespace blitz_ext

#endif //AT_BZUTIL_H__
