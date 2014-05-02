//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_BZARRAY_ITERATOR_H__
#define AT_BZARRAY_ITERATOR_H__

//As of version of Blitz (09/99) there is no support for
//iterator_traits<ArrayIterator>. This header adds that
//using partial specialization.

#if 0

#include <blitz/array.h>

#ifdef _STLP_STD
namespace _STLP_STD
#else
namespace std
#endif
{

template <class T, int N>
struct iterator_traits<blitz::ConstArrayIterator<T,N> > {
  //Note: iterator_category must be random_access_iterator,
  //but iterator implementation does not yet provide
  //necessary operators (like -()(iter,iter) etc).
  typedef forward_iterator_tag      iterator_category;
  typedef  T                              value_type;
  typedef int                             difference_type;
  typedef const T*                        pointer;
  typedef const T&                        reference;
};

template <class T, int N>
struct iterator_traits<blitz::ArrayIterator<T,N> > {
  //see note above for iterator_category choice.
  typedef forward_iterator_tag      iterator_category;
  typedef  T                              value_type;
  typedef int                             difference_type;
  typedef T*                              pointer;
  typedef T&                              reference;
};

} //namespace std;

#endif

#endif //AT_BZARRAY_ITERATOR_H__
