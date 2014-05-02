//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/Math/pfactors.hpp"

#include "PRODDL/Math/fft.hpp"

#include "PRODDL/Testing/exception.hpp"

#include "PRODDL/Common/debug.hpp"

#include "gtest/gtest.h"

namespace {

int powInt(int a, int exp) {
  int res = 1;
  for(int i = 0; i < exp; i++) res *= a;
  return res;
}

int factToNumber(const int *factors, const int *exp, int n) {
  int res = 1;
  for( int i = 0; i < n; i++ ) res *= powInt(factors[i], exp[i]);
  return res;
}

// set numbers[i] to 1 if i can be factored into primes[], set other numbers[] to 0.
// complexity: nNumbers^2

void seeve(int numbers[], const int nNumbers, const int primes[], const int nPrimes) {
  // reset numbers[] to 0
  for(int i=0; i < nNumbers; i++)
    numbers[i] = 0;
  // seed numbers[] with primes[]
  for(int i = 0; i < nPrimes; i++) {
    int prime = primes[i];
    if(prime < nNumbers) {
      numbers[prime] = 1;
    }
  }
  // set all multiplies of already found elements in numbers[] to 1
  for(int i = 2; i < nNumbers; i++) {
    if(numbers[i]) {
      for(int j = 2; j <= i; j++) {
	if(numbers[j]) {
	  int product = i*j;
	  if(product < nNumbers)
	    numbers[product] = 1;
	}
      }
    }
  }
}

int findFirstFactored(int numbers[],int nNumbers,int i_start) {
  for( ; i_start < nNumbers; i_start++) {
    if(numbers[i_start])
      return i_start;
  }
  return i_start;
}

} // namespace


TEST(PfactorsTest, All) {

  using namespace PRODDL::Math;

  const int primes[] = { 2, 3, 5, 7, 11, 13 };

  const int nPrimes = sizeof(primes)/sizeof(primes[0]);

  const int ind11 = 4, ind13 = 5;

  const int nNumbers = 1024;

  int numbers[nNumbers];

  seeve(numbers,nNumbers,primes,nPrimes);

  int exponents[nPrimes];

  for(int i = 0; i < nNumbers; i++) {
    int nextFactored = findFirstFactored(numbers,nNumbers,i);
    int nFactoredFound = 0;
    bool is_found = findFirstPrimeFactoring(i,nNumbers,primes,primes+nPrimes,exponents,nFactoredFound);
    //DEBUG:
    ATOUTVAR(nextFactored); ATOUTVAR(is_found); ATOUTVAR(nFactoredFound); 
    ATOUTRANGE(exponents,exponents+nPrimes,int); ATOUTENDL();
    if(nextFactored < nNumbers) {
      if( ! is_found ) {
	throw PRODDL::Testing::test_error("Did not find expected factorization.");
      }
      else {
	if(nFactoredFound != nextFactored)
	  throw PRODDL::Testing::test_error("Found factored number is not equal to the expected one.");
	else {
	  int nFromExponents = factToNumber(primes, exponents, nPrimes);
	  if( nFromExponents != nFactoredFound ) {
	    throw PRODDL::Testing::test_error("Exponents do not correspond to the factored number.");
	  }
	  if(exponents[ind11] + exponents[ind13] <= 1) {
	    ATOUTVAR(i); ATOUTVAR(nFactoredFound);
	    int iFFTW_Size = 0;
	    if( FFTW_Size::findBestSize(i,iFFTW_Size)) {
	      ATOUTVAR(iFFTW_Size);
	      if( nFactoredFound != iFFTW_Size ) {
		throw PRODDL::Testing::test_error("Found size for FFTW is not equal to the expected one.");
	      }
	    }
	    else {
	      throw PRODDL::Testing::test_error("FFTW best size could not be found.");
	    }
	    ATOUTENDL();
	  }
	}
      }
    }
    else if(is_found) {
      throw PRODDL::Testing::test_error("Unexpected factorization found.");
    }
  }

}
  
