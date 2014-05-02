//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_BZATIO_H__
#define AT_BZATIO_H__

//i/o extensions to Blitz

#include <blitz/array.h>

namespace blitz_ext
{

  using namespace blitz;

  // Print array making each row to have number of elements equal to the size of
  // the last (fastest) dimension of array. So, 2-D array size [2 x 3] will look like:
  // a11 a12 a13
  // a21 a22 a23
  // (contast to Blitz which prints max seven elements on a row).
  // Dimensions are separated by the number of newline characters equal to
  // the rank of dimension.
  // Array is printed in storage order (like in Blitz operator<< ).

  template<class T_numtype, int N_rank>
    std::ostream& PrintAsMatrix(std::ostream& os, const Array<T_numtype,N_rank>& x)
    {
      for (int i=0; i < N_rank; ++i)
	{
	  os << x.extent(i);
	  if (i != N_rank - 1)
	    os << " x ";
	}
      
      os << std::endl << "[ ";
      
      typename Array<T_numtype, N_rank>::const_iterator iter = x.begin(), end = x.end();

      int p = 0;
      while (iter != end) {
	
	// See if we need a linefeed and how many

	TinyVector<int,N_rank> offset_ind;
        offset_ind = iter.position() - x.base();

	if(p > 0)
	  {
	    //if offset_ind is, say, [ 2, 3, 0, 0 ] and array is stored in C-order,
	    //the following loop will print two newline characters

	    for(int i=0; i < N_rank && (offset_ind(x.ordering(i)) == 0); i++) 
	      os << std::endl << "  ";
	  }

	// print the data element

	os << std::setw(9) << (*iter) << " ";

	++p;
	++iter;
      }
      
      os << "]" << std::endl;
      return os;
    }  
  
  
} //namespace blitz_ext

#endif //AT_BZATIO_H__
