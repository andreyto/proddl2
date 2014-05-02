//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_MATH_FFT_H_
#define AT_MATH_FFT_H_

// Functions related to FFT: finding the best grid size etc

#include <iterator> // for std::advance()

#include "PRODDL/Math/pfactors.hpp"

#include "PRODDL/Math/math_except.hpp"

#include "PRODDL/Common/common_types.hpp"

#include "PRODDL/Common/nd_index_iter.hpp"

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/logger.hpp"

namespace PRODDL { namespace Math {


  // Class that defines static methods to compute "best" (fastest) sizes for FFTW based on
  // factoring into primes preffered by FFTW.

  class FFTW_Size {

  public:

  protected:

    // The number of prime factors

    enum { _nPrimes = 6 };

    // The prime factors that provide best sizes for this FFTW.

    static const int _primes[_nPrimes]; 


  public:


    // Return the number of prime factors used by this FFT

    static inline int getNumberOfPrimes() {
      return _nPrimes;
    }


    // Return the array of prime factors used by this FFT (size of array is provided by 'getNumberOfPrimes()')

    static inline const int* getPrimes() {
      return _primes;
    }


    // Find first positive integer number located in range ['startN','endN') that can be factored into
    // a subset of primes most efficiently handled by FFTW Fast Fourier Transform library.
    // This function does the same as PRODDL::Math::findFirstPrimeFactoring, but uses the first 'nPrimes' from a fixed set of
    // prime factors [2,3,5,7,11,13] with the additional condition for the exponents in factoring 
    // that if N = 2^i*...*11^j*13^l then j+l <= 1. Documentation for FFTW library v.2.1.3 states that
    // array sizes with such prime factoring are handled most efficiently (although they say that the additional source code
    // for FFTW can be generated to provide custom processing for other prime factors).
    // Requirements: 'nPrimes' must be <= 'getNumberOfPrimes()'. If 'nPrimes' is omitted,
    // the default value of 'getNumberOfPrimes()' will be used (that is probably how you will call this function most of the time).
    // Results: write exponents of prime factors
    // into the forward iterator range ['resultExponent','resultExponent' + nPrimes).
    // 'resultExponent' is required to be forward iterator model because we repeatedly use it as a temporary work buffer
    // inside this function. The found factored number is returned in 'resultN' argument. Function returns
    // true if the factoring is successfully found, and false otherwise. In the later case the content of the
    // output arguments is undefined.
  
    template<class ForwardIterExponent> 
    static bool findFirstPrimeFactoring(int startN, int endN,
					ForwardIterExponent resultExponent,
					int& resultN,
					int nPrimes = _nPrimes) throw(fft_size_error)
    {
    
      if( nPrimes > _nPrimes || nPrimes < 0 ) {
	throw fft_size_error("FFTW_Size::findFirstPrimeFactoring(): "
			     "nPrimes argument must be less than getNumberOfPrimes() and not less than zero.");
      }

      // indexes of primes in '_primes' array, which have special condition imposed on them
      const int ind11 = 4, ind13 = 5;

      // create iterators pointing to resultExponent elements for 11 and 13 - we will then dereference
      // them repeatedly within a loop
      ForwardIterExponent iterExp11 = resultExponent;
      std::advance(iterExp11,ind11);
      ForwardIterExponent iterExp13 = resultExponent;
      std::advance(iterExp13,ind13);

      while(1) {
	if( Math::findFirstPrimeFactoring(startN,endN,_primes,_primes+nPrimes,resultExponent,resultN) ) {
	  if(*iterExp11 + *iterExp13 <= 1) {
	    return true;
	  }
	  startN = resultN + 1;
	}
	else {
	  return false;
	}
      }
    }

    // Calls findFirstPrimeFactoring() function with 'endN' argument set to the maximum value of
    // a given int as defined by std::numeric_limits.

    template<class ForwardIterExponent> 
    static bool findFirstPrimeFactoring(int startN,
					ForwardIterExponent resultExponent,
					int& resultN,
					int nPrimes = _nPrimes) 
    {
      const int maxN = std::numeric_limits<int>::max();
      return findFirstPrimeFactoring(startN, maxN,
				     resultExponent,
				     resultN,
				     nPrimes);
    }

    // Set 'resultN' to the first size for FFT such that 'resultN' >= 'startN'.
    // Return true if found, false - otherwise.
    // If 'nPrimes' is set, only first 'nPrimes' elements of 'getPrimes()' array will be used to
    // search for factoring, which determines the best FFT size.
    // This function calls 'findFirstPrimeFactoring()' but ignores found exponents.
    // Note: from raw files with FFTW benchmark results, it appears that almost always
    // the smaller size is faster, regardless of the prime factors (e.g. 60 is faster than 64).

    static bool findBestSize(int startN,
			     int& resultN,
			     int nPrimes = _nPrimes) 
    {
      int exponents[_nPrimes];
      const int maxN = std::numeric_limits<int>::max(); 
      return findFirstPrimeFactoring(startN, maxN,
				     exponents,
				     resultN,
				     nPrimes);
    }


    // multidimensional version

    template<int n_dim>
    static
    typename common_types::point_type<bool,n_dim>::Type
    findBestSize(const typename common_types::point_type<int,n_dim>::Type& startN,
		 typename common_types::point_type<int,n_dim>::Type& resultN,
		 int nPrimes = _nPrimes)
    {
      typename common_types::point_type<bool,n_dim>::Type ret;
      for(int i = 0; i < n_dim; i++) {
	ret(i) = findBestSize(startN(i),resultN(i),nPrimes);
      }
      return ret;
    }

  }; // class FFTW_Size


  // class to handle wrapped-around periodic results produced by application
  // of FFT

  template<typename T_num,int n_dim>
  class WrappedIndex {

  public:

    typedef typename common_types::point_type<int,n_dim>::Type IntPoint;

    typedef typename common_types::num_multiarray_type<T_num,n_dim>::Type Array;

    typedef index_mover<n_dim> ind_mover;

    WrappedIndex() :
      size(0), neg_half_size(0), lbound(0)
    {
    }

    WrappedIndex(const IntPoint& _size) {
      // Only zero-bazed arrays are supported.
      lbound = 0;
      size = _size;
      neg_half_size = -1 * (size/2);
    }

    IntPoint unwrap(const IntPoint& x) const {
      IntPoint y = -1 * x;
      for(int i=0; i < n_dim; i++) {
	if(y(i) < neg_half_size(i))
	  y(i) += size(i);
      }
      return y;
    }


    //TODO: unwrapping arrays can be made for efficient by
    //moving entire quadrants of data instead of iterating
    //over the index


    // Copy data from 'in' array into 'out' array in unwrapped
    // order. No matter what the shape of the arrays is,
    // only data in [0,size) range will be copied.
    // The unwrapped data is centered at the middle of the
    // output array, so that index (size/2) holds frequency zero.

    void
    unwrap(const Array& in, Array& out) {

      dbg::trace t1(DBG_HERE);

      // Doing it on the same data as input and output will
      // give wrong results.
      // The method we use below to check that data is different
      // will not work if arrays only partly overlap.
      // The assertion inside the inner loop makes the universal
      // check that no two elements of 'in' and 'out' share
      // the same memory location.

      dbg::assertion(dbg::error, DBG_ASSERTION(in.dataFirst() != out.dataFirst()));

      // ... Loop over all indices ...

      IntPoint ind;

      for( ind_mover::before_first(ind,lbound); 
	   ind_mover::next(ind,lbound,size);  ) {
	
	IntPoint indUnwr = unwrap(ind) - neg_half_size;

	out(indUnwr) = in(ind);

	dbg::assertion(dbg::error, DBG_ASSERTION( &(const_cast<Array&>(in)(ind)) != &out(ind)));

      }
      
    }

    IntPoint before_first() const {
      IntPoint ind;
      ind_mover::before_first(ind,lbound);
      return ind;
    }

    bool next(IntPoint& ind) const {

      return ind_mover::next(ind,lbound,size);

    }

  protected:

    IntPoint size;
    IntPoint neg_half_size;
    IntPoint lbound;

  }; // class WrappedIndex<>
  
} } // namespace PRODDL { namespace Math {

#include "PRODDL/Math/fft.cc"

#endif //  AT_MATH_FFT_H_ 
