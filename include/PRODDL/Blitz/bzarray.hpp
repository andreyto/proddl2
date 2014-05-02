//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_BLITZ_BZARRAY_H__
#define PRODDL_BLITZ_BZARRAY_H__

#include <blitz/array.h>

// Some utility classes and functions for Blitz Array

namespace PRODDL {


  // Class the describes the structure of Blitz Array (i.e. lower bounds, extents, storage, strides).
  // It is able to load this data from the existing Array and create a new Array with
  // a specified structure (when a new array is created, strides are computed automatically
  // from other parameters, because a new array is always contiguous in memory).
  // Rational: in parallel programs, we need to pass first the structure of the array in
  // order to create a receiving array on the other side of the communication.

template<int n_dim>
class ArrayStruct {
public:

  typedef ArrayStruct<n_dim> Self;

  typedef blitz::TinyVector<int,n_dim> Dimensions;

  typedef blitz::GeneralArrayStorage<n_dim> Storage;

  ArrayStruct() {}

  template<typename T_num>
  ArrayStruct(const blitz::Array<T_num,n_dim>& arr) {

    init(arr);

  }
  
  template<typename T_num>
  void
  init(const blitz::Array<T_num,n_dim>& arr) {

    extents = arr.extent();
    strides = arr.stride();

   for(int i = 0; i < n_dim; i++) {
     storage.ascendingFlag()(i) = arr.isRankStoredAscending(i);
   }

   storage.base() = arr.base();
   storage.ordering() = arr.ordering();
  }

  template<typename T_num> 
  blitz::Array<T_num,n_dim>
  createArray(T_num) const {

    return blitz::Array<T_num,n_dim>(storage.base(),extents,storage);

  }

  //TODO: Add access methods

  bool operator== (const Self& other) {

    return blitz::all(extents == other.extents) &&
      blitz::all(strides == other.strides) &&
      blitz::all(storage.base() == other.storage.base()) &&
      blitz::all(storage.ordering() == other.storage.ordering()) &&
      blitz::all(ascendingFlag() == other.storage.ascendingFlag());

  }

  template<typename T_num>
  bool isSameStruct(const blitz::Array<T_num,n_dim>& arr)
  {

    bool one =  blitz::all(extents == arr.extent()) &&
      blitz::all(strides == arr.stride()) &&
      blitz::all(storage.base() == arr.base()) &&
      blitz::all(storage.ordering() == arr.ordering());
    bool two = true;
    for(int i = 0; i < n_dim; i++) {
      if( storage.ascendingFlag()(i) != arr.isRankStoredAscending(i) )
	{
	  two = false;
	}
    }
    return one && two;
  }

    
protected:

  Dimensions extents;
  Dimensions strides;
  Storage storage;

};


  // specialization of blitz::GeneralArrayStorage for
  // row-major (c-style) arrays with possibly non-zero base

  template<int N_rank>
  class RowMajorArrayAnyBase : public blitz::GeneralArrayStorage<N_rank> {
  public:
    RowMajorArrayAnyBase(const blitz::TinyVector<int,N_rank>& _base)
      : blitz::GeneralArrayStorage<N_rank>()
    {
      base_ = _base;
    }
  };


} // namespace PRODDL

#endif // PRODDL_BLITZ_BZARRAY_H__
