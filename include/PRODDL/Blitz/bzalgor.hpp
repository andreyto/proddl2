//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_BZALGOR_H__
#define AT_BZALGOR_H__

//some extensions to Blitz

#include <blitz/tinyvec2.h>

#include "PRODDL/Common/math.hpp"

namespace blitz_ext
{

namespace internal
  {

template<class T,int N,int I> struct recursion_tinyvec;

template<class T, int N> struct recursion_tinyvec<T,N,N>
{
  inline static void
  bracket(const blitz::TinyVector<T,N>& x,
	      blitz::TinyVector<T,N>& x_min, 
	      blitz::TinyVector<T,N>& x_max) { }
  inline static void
  crop(blitz::TinyVector<T,N>& x,
       const blitz::TinyVector<T,N>& x_min, 
       const blitz::TinyVector<T,N>& x_max,
       bool& cropped) { }
  inline static bool
  in_range(const blitz::TinyVector<T,N>& x,
       const blitz::TinyVector<T,N>& x_min, 
       const blitz::TinyVector<T,N>& x_max) { return true; }
  inline static bool
  between(const blitz::TinyVector<T,N>& x,
	  const blitz::TinyVector<T,N>& x_min, 
	  const blitz::TinyVector<T,N>& x_max) { return true; }
  template<class Predicate2> inline static bool
  is_each(const blitz::TinyVector<T,N>& x1,
	  const blitz::TinyVector<T,N>& x2, 
	  Predicate2 pred) { return true; }
  template<class Predicate2> inline static bool
  is_any(const blitz::TinyVector<T,N>& x1,
	 const blitz::TinyVector<T,N>& x2, 
	 Predicate2 pred) { return false; }
  inline static void
  min_each(const blitz::TinyVector<T,N>& x1, 
	   const blitz::TinyVector<T,N>& x2,
	   blitz::TinyVector<T,N>& x_res) { }
  inline static void
  max_each(const blitz::TinyVector<T,N>& x1, 
	   const blitz::TinyVector<T,N>& x2,
	   blitz::TinyVector<T,N>& x_res) { }
};

template<class T,int N,int I> struct recursion_tinyvec
{
  //find min and max values for each component of x:
  inline static void
    bracket(const blitz::TinyVector<T,N>& x,
	    blitz::TinyVector<T,N>& x_min, 
	    blitz::TinyVector<T,N>& x_max) {
    //Note: we use the 'else' between
    //two 'if' statements below, 
    //so bracket(1,100,-100) will give (1,-100).
    //Such behavior requires initialization of
    //x_min and x_max with the first element when
    //bracket() is iterated over the range.
    if(x[I] < x_min[I]) x_min[I] = x[I];
    else if(x[I] > x_max[I]) x_max[I] = x[I];
    recursion_tinyvec<T,N,I+1>::bracket(x,x_min,x_max);
  }
  //if (x < x_min) { x = x_min; cropped = true; },
  //same for x_max. Requirements: x_min <= x_max.
  inline static void
    crop(blitz::TinyVector<T,N>& x,
	 const blitz::TinyVector<T,N>& x_min, 
	 const blitz::TinyVector<T,N>& x_max,
	 bool& cropped) {
    if(x[I] < x_min[I]) { x[I] = x_min[I]; cropped = true; }
    else if(x[I] > x_max[I]) { x[I] = x_max[I]; cropped = true; }
    recursion_tinyvec<T,N,I+1>::crop(x,x_min,x_max,cropped);
  }
  // Return true if x is in [x_min, x_max) in all dimensions.
  // Requirements: x_min <= x_max.
  inline static bool
    in_range(const blitz::TinyVector<T,N>& x,
	     const blitz::TinyVector<T,N>& x_min, 
	     const blitz::TinyVector<T,N>& x_max) {
    if(x[I] < x_min[I]) { return false; }
    else if(x[I] >= x_max[I]) { return false; }
    else return recursion_tinyvec<T,N,I+1>::in_range(x,x_min,x_max);
  }
  //same as "in_range" but right boundary inclusive
  inline static bool
    between(const blitz::TinyVector<T,N>& x,
	    const blitz::TinyVector<T,N>& x_min, 
	    const blitz::TinyVector<T,N>& x_max) {
    if(x[I] < x_min[I]) { return false; }
    else if(x[I] > x_max[I]) { return false; }
    else return recursion_tinyvec<T,N,I+1>::between(x,x_min,x_max);
  }
  //apply binary predicate pred to each pair of corresponding
  //elements from x1 and x2 and return true if predicate yeld
  //true for all pairs
  template<class Predicate2> inline static bool
  is_each(const blitz::TinyVector<T,N>& x1,
	  const blitz::TinyVector<T,N>& x2, 
	  Predicate2 pred) { 
    if( ! pred(x1[I],x2[I])) { return false; }
    else return recursion_tinyvec<T,N,I+1>::is_each(x1,x2,pred);
  }
  //apply binary predicate pred to each pair of corresponding
  //elements from x1 and x2 and return true if predicate yeld
  //true for any pair
  template<class Predicate2> inline static bool
  is_any(const blitz::TinyVector<T,N>& x1,
	 const blitz::TinyVector<T,N>& x2, 
	 Predicate2 pred) { 
    if( pred(x1[I],x2[I])) { return true; }
    else return recursion_tinyvec<T,N,I+1>::is_any(x1,x2,pred);
  }
  //for each I, set x_res[I] = min(x1[I],x2[I])
  inline static void
  min_each(const blitz::TinyVector<T,N>& x1, 
	   const blitz::TinyVector<T,N>& x2,
	   blitz::TinyVector<T,N>& x_res) {
    if(x1[I] < x2[I]) x_res[I] = x1[I];
    else x_res[I] = x2[I];
    recursion_tinyvec<T,N,I+1>::min_each(x1,x2,x_res);
  }
  //for each I, set x_res[I] = max(x1[I],x2[I])
  inline static void
  max_each(const blitz::TinyVector<T,N>& x1, 
	   const blitz::TinyVector<T,N>& x2,
	   blitz::TinyVector<T,N>& x_res) {
    if(x2[I] < x1[I]) x_res[I] = x1[I];
    else x_res[I] = x2[I];
    recursion_tinyvec<T,N,I+1>::max_each(x1,x2,x_res);
  }
};

  } //namespace internal

template<class T,int N> 
void bracket(const blitz::TinyVector<T,N>& x,
	blitz::TinyVector<T,N>& x_min, blitz::TinyVector<T,N>& x_max)
{
  internal::recursion_tinyvec<T,N,0>::bracket(x,x_min,x_max);
}

template<class T,int N> 
void crop(blitz::TinyVector<T,N>& x,
	  const blitz::TinyVector<T,N>& x_min, 
	  const blitz::TinyVector<T,N>& x_max,
	  bool& cropped)
{
  cropped = false;
  internal::recursion_tinyvec<T,N,0>::crop(x,x_min,x_max,cropped);
}

template<class T,int N> 
bool in_range(const blitz::TinyVector<T,N>& x,
	      const blitz::TinyVector<T,N>& x_min, 
	      const blitz::TinyVector<T,N>& x_max)
{
  return internal::recursion_tinyvec<T,N,0>::in_range(x,x_min,x_max);
}

template<class T,int N> 
bool between(const blitz::TinyVector<T,N>& x,
	     const blitz::TinyVector<T,N>& x_min, 
	     const blitz::TinyVector<T,N>& x_max)
{
  return internal::recursion_tinyvec<T,N,0>::between(x,x_min,x_max);
}

template<class T,int N, class Predicate2> 
bool is_each(const blitz::TinyVector<T,N>& x1,
	     const blitz::TinyVector<T,N>& x2, 
	     Predicate2 pred)
{
  return internal::recursion_tinyvec<T,N,0>::is_each(x1,x2,pred);
}

template<class T,int N, class Predicate2> 
bool is_any(const blitz::TinyVector<T,N>& x1,
	    const blitz::TinyVector<T,N>& x2, 
	    Predicate2 pred)
{
  return internal::recursion_tinyvec<T,N,0>::is_any(x1,x2,pred);
}

template<class T,int N> 
void  min_each(const blitz::TinyVector<T,N>& x1, 
	       const blitz::TinyVector<T,N>& x2,
	       blitz::TinyVector<T,N>& x_res)
{
  return internal::recursion_tinyvec<T,N,0>::min_each(x1,x2,x_res);  
}

template<class T,int N> 
void  max_each(const blitz::TinyVector<T,N>& x1, 
	       const blitz::TinyVector<T,N>& x2,
	       blitz::TinyVector<T,N>& x_res)
{
  return internal::recursion_tinyvec<T,N,0>::max_each(x1,x2,x_res);  
}

//function object adaptors:

struct bracketer
{
  template<class T,int N> void operator() (const blitz::TinyVector<T,N>& x,
		   blitz::TinyVector<T,N>& x_min, blitz::TinyVector<T,N>& x_max) {
    internal::recursion_tinyvec<T,N,0>::bracket(x,x_min,x_max);
  }
};

struct cropper
{
  template<class T,int N> void operator() (blitz::TinyVector<T,N>& x,
					   const blitz::TinyVector<T,N>& x_min, 
					   const blitz::TinyVector<T,N>& x_max,
					   bool& cropped) {
    internal::recursion_tinyvec<T,N,0>::crop(x,x_min,x_max,cropped);
  }
};

} //namespace blitz_ext

#endif //AT_BZALGOR_H__
