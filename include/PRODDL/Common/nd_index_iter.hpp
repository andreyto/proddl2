//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef  AT_ND_INDEX_ITER_H_
# define  AT_ND_INDEX_ITER_H_

// Generic routine to iterate multidimensional index
// Example:
//   int array[3][3][4];
//   int lbound[3] = {0,0,0};
//   int ubound[3] = {3,3,4};
//   int ind[3];
//   typedef index_mover<3> ind_mover;
//   for(ind_mover::first(ind,lbound); ind_mover::is_good(ind,ubound); ind_mover::next(ind,lbound,ubound)) {
//      array[ind[0]][ind[1]][ind[2]] = ...;
//   }
//
//  ...
// Slightly more efficient use would be (easier for the compiler to optimize out the redundant 'if' construct):
// 
//   for(ind_mover::before_first(ind,lbound); ind_mover::next(ind,lbound,ubound); ) {
//      array[ind[0]][ind[1]][ind[2]] = ...;
//   }
//
// The idea is to have dimension-independent notation to iterate through indices of multidimensional array -
// useful for array types that have dimension-independent notation for indexing (not like simple example above
// with c-array). This is a primitive implementation - for the real stuff see, for instance, Blitz lib: <blitz/array/eval.cc>.
// This is tuned currently for performance on row-major arrays - the last dimension of index will be changed the fastest,
// although it will still work correctly on other orderings (e.g. Fortran) albeit inefficiently.
// Requirements: Ind must have indexing operator[] defined. n_dim must be > 0.
// Performance: is_good() makes one comparison; next() makes one comparison 'most of the time' (in the fastest
// changing dimension). Thus, w/o compiler optimization it is typically one comparison more than normal multi-level loop syntaxis.
// Postconditions: after loop finished, 'ind' points to the first element after the allowed domain.
//
// ATTENTION: methods here expect half-open range, like STL iterators. If using Blitz Array ubound() (which is a closed range), 
// use next(ind,arr.lbound(),arr.ubound()+1).
//
// TODO: create partial specialization for 1-D case; introduce orderings other than row-major.

namespace PRODDL {

template<int n_dim> class index_mover {
 
public:

  enum { max_dim = n_dim - 1 };

  // advance index 'ind' to the next position;
  // return 'true' if the new index is within bounds, false - otherwise
  template<class Ind> inline static bool next(Ind& ind, 
					      const Ind& lbound,
					      const Ind& ubound) {
    
    // we hit this most of the time
    if( ++ind[max_dim] < ubound[max_dim] ) return true;
    
    // reset previous dimensions to lbound, increase next dimension and check again
    for(int i_dim = max_dim-1; i_dim >= 0; i_dim--) {
      ind[i_dim+1] = lbound[i_dim+1];
      if( ++ind[i_dim] < ubound[i_dim] ) return true;
    }
    // if we got here, it means 'ind[0] >= ubound[0]' - time to stop
    return false;
  }
  
  template<class Ind> inline static void first(Ind& ind, 
					       const Ind& lbound) {
    for ( int i = 0; i < n_dim; i++) ind[i] = lbound[i]; 
  }

  // set 'ind' to such value, that the next call to 'next()'
  // will set 'ind' to the first valid index value.
  // use this function to initialize 'ind' in loops like this
  // for(before_first(ind,lbound); next(ind,lbound,ubound); ) {...}
  template<class Ind> inline static void before_first(Ind& ind, 
						      const Ind& lbound) {
    if(n_dim > 0) ind[max_dim] = lbound[max_dim] - 1;
    for ( int i = 0; i < max_dim; i++) ind[i] = lbound[i];
  }
  
  template<class Ind> inline static bool is_good(const Ind& ind, 
						 const Ind& ubound) {
    return ( ind[0] < ubound[0] ); 
  }

};

} // namespace PRODDL

#endif // AT_ND_INDEX_ITER_H_
