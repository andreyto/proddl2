//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_MATH_PFACTOR_H_
#define AT_MATH_PFACTOR_H_

// Functions related to prime factoring

#include <limits>

namespace PRODDL { namespace Math {

  // Find first positive integer number located in range ['startN','endN') that can be factored into
  // a set of primes from iterator range ['firstPrime', 'lastPrime'). Write exponents of prime factors
  // into the forward iterator range ['resultExponent','resultExponent' + ('lastPrime' - 'firstPrime')).
  // 'resultExponent' is required to be forward iterator model because we repeatedly use it as a temporary work buffer
  // inside this function. The found factored number is returned in 'resultN' argument. Function returns
  // true if the factoring is successfully found, and false otherwise. In the later case the content of the
  // output arguments is undefined.
  // This function uses brute force approach of repeatedly dividing by prime factors, so it is suited for relatively 
  // small problems. There are more sofisticated algorithms in Netlib. Originally, this function was designed to
  // find the smallest grid size for FFT.

  template<class ForwardIterExponent, class ForwardIterPrime, typename IntegerType> 
  bool findFirstPrimeFactoring(IntegerType startN, IntegerType endN,
			       ForwardIterPrime firstPrime, ForwardIterPrime lastPrime,
			       ForwardIterExponent resultExponent,
			       IntegerType& resultN)
  {
    
    for(IntegerType N = ( startN <= 1 ? 2 : startN ); N < endN; ++N) {
      ForwardIterExponent pExponent = resultExponent;
      IntegerType n = N;
      for( ForwardIterPrime pPrime = firstPrime; pPrime != lastPrime; ++pPrime, ++pExponent ) {
	IntegerType prime = *pPrime;
	IntegerType exponent = 0;
	while(n % prime == 0) {
	  n /= prime;
	  ++exponent;
	}
	*pExponent = exponent;
	if( n == 1 ) { // found factoring, now fill the remaining exponents with zeros and return
	  for( ++pPrime, ++pExponent ; pPrime != lastPrime; ++pPrime, ++pExponent ) {
	    *pExponent = IntegerType(0);
	  }
	  resultN = N;
	  return true;
	}
      } // for( ForwardIterPrime pPrime
    } // for(IntegerType N
    return false;
  }

  // Calls findFirstPrimeFactoring() function with 'endN' argument set to the maximum value of
  // a given IntegerType as defined by std::numeric_limits.

  template<class ForwardIterExponent, class ForwardIterPrime, typename IntegerType> 
  bool findFirstPrimeFactoring(IntegerType startN,
			       ForwardIterPrime firstPrime, ForwardIterPrime lastPrime,
			       ForwardIterExponent resultExponent,
			       IntegerType& resultN) 
  {
    const IntegerType maxN = std::numeric_limits<IntegerType>::max();
    return findFirstPrimeFactoring(startN, maxN,
				   firstPrime, lastPrime,
				   resultExponent,
				   resultN);
  }
  
} } // namespace PRODDL { namespace Math {

#endif //  AT_MATH_PFACTOR_H_ 
