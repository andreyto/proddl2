//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_COMMON_MATH_H__
# define AT_COMMON_MATH_H__

// This included header must be standard conforming
#include <cmath>

namespace PRODDL { namespace Math {

  // Small powers. Adapted from Blitz.

  template<typename T>  inline T pow2(T x) { return x*x; }
  template<typename T>  inline T pow3(T x) { return x*x*x; }
  template<typename T>  inline T pow4(T x) { T t1 = x*x; return t1*t1; }
  template<typename T>  inline T pow5(T x) { T t1 = x*x; return t1*t1*x; }
  template<typename T>  inline T pow6(T x) { T t1 = x*x*x; return t1*t1; }
  template<typename T>  inline T pow7(T x) { T t1 = x*x; return t1*t1*t1*x; }
  template<typename T>  inline T pow8(T x) { T t1 = x*x, t2=t1*t1; return t2*t2; }
  template<typename T>  inline T pow12(T x) { return pow2(pow6(x)); }

  template<typename T>  inline T sign(T x, T y) { return y >= 0 ? std::abs(x) : -std::abs(x); }

  //TODO: Commented out static initialization code is broken - it depends on the order of the instantiation of static members.
  //Should be fixed using Miers Singleton concept
  

  template<typename T_num>
  struct Constants {
    typedef T_num num_type;

    typedef Constants Self;

    static T_num PI;
    static T_num TWO_PI;
    static T_num HALF_PI;
    static T_num DegToRad;

    const static T_num Angstrom = 1e-10;
    const static T_num Nano = 1e-9;
    const static T_num AngstromToNm = 0.1; // angstrom in nanometers

    static Constants statSelf;

    Constants() {

      Self::PI = T_num(4.0)*std::atan(T_num(1.0));
      Self::TWO_PI = T_num(2.0)*Self::PI;
      Self::HALF_PI = T_num(0.5)*Self::PI;
      Self::DegToRad = Self::PI/T_num(180.);

    }

  };

  template<typename T_num> T_num Constants<T_num>::PI;
  template<typename T_num> T_num Constants<T_num>::TWO_PI;
  template<typename T_num> T_num Constants<T_num>::HALF_PI;
  template<typename T_num> T_num Constants<T_num>::DegToRad;
  template<typename T_num> Constants<T_num> Constants<T_num>::statSelf;


//   template<typename T_num> const T_num Constants<T_num>::PI = T_num(4.0)*std::atan(T_num(1.0));
//   template<typename T_num> const T_num Constants<T_num>::TWO_PI = T_num(2.0)*Constants<T_num>::PI;
//   template<typename T_num> const T_num Constants<T_num>::HALF_PI = T_num(0.5)*Constants<T_num>::PI;


//   template<typename T_num> const T_num Constants<T_num>::DegToRad = Constants<T_num>::PI/T_num(180.);

  // This function tests whether or not x and y of an 
  // integer or real type are equal subject to the given relative and absolute tolerances. 
  // The formula used is:
  //
  // | x - y | < atol + rtol * | y |
  // The definition is taken from allclose() function in Python Numarray

  template<typename T_num1, typename T_num2, typename T_num3> 
  bool
  are_close (T_num1 x, T_num2 y, T_num3 rtol, T_num3 atol) 
  {
    return std::abs( x - y ) < atol + rtol * std::abs( y );
  }

  template<typename T_num1, typename T_num2> 
  bool
  are_close (T_num1 x, T_num2 y) 
  {
    return std::abs( x - y ) < 1.e-8 + 1.e-5 * std::abs( y );
  }


  // Synonym for are_close()

  template<typename T_num1, typename T_num2, typename T_num3> 
  bool are_equal_delta(T_num1 x, T_num2 y, T_num3 delta) {
    return are_close(x,y,delta,delta);
  }

} } // namespace PRODDL::Math

#endif // AT_COMMON_MATH_H__
