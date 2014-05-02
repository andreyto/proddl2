//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_BZTINYVEC_H_
#define AT_BZTINYVEC_H_

//Some extensions to blitz::TinyVector and blitz::TinyMatrix

#include <blitz/array.h>
#include <blitz/tinymat2.h>
#include <blitz/tinymat2.cc>

#include "PRODDL/Blitz/bzalgor.hpp"

#include "PRODDL/Common/math.hpp"

#include <functional>

//specialization of equal_to functor for TinyVector that returns
//a scalar bool by AND-combining results of pairwise comparison.
//Note: the operator== defined by Blitz returns an array of bool.
//Arguable: inject it into std namespace for it to work in all cases

#ifdef _STLP_STD
namespace _STLP_STD {
#else
namespace std {
#endif

	template <class T_value, int N_elem>
	struct equal_to<blitz::TinyVector<T_value,N_elem> > : 
		public binary_function<blitz::TinyVector<T_value,N_elem>,blitz::TinyVector<T_value,N_elem>,bool> 
	{
		bool operator()(const blitz::TinyVector<T_value,N_elem>& x1, 
			const blitz::TinyVector<T_value,N_elem>& x2) const 
		{ return blitz_ext::is_each(x1,x2,std::equal_to<T_value>()); }
	};

} //namespace std



namespace blitz_ext
{
	/// matrix-vector product for tiny matrix and tiny vector
	// Somehow this function signature has been removed from Blitz v0.10
	template<typename T, int N_rows, int N_columns>
	blitz::TinyVector<T,N_rows>
		product(
		const blitz::TinyMatrix<T,N_rows,N_columns>& m,
		const blitz::TinyVector<T,N_columns>& x
		)
	{
		blitz::TinyVector<T,N_rows> y;
		for(int i=0; i < N_rows; ++i) {
			y(i) = 0;
			for(int j=0; j < N_columns; ++j)
				y(i) += m(i,j)*x(j);
		}
		return y;
	}

	/// matrix-matrix product for tiny matrix and tiny matrix
	// Somehow this function signature has been removed from Blitz v0.10
	template<typename T, int N_rows, int N_columns, int M_columns>
	blitz::TinyMatrix<T,N_rows,M_columns>
		product(
		const blitz::TinyMatrix<T,N_rows,N_columns>& x,
		const blitz::TinyMatrix<T,N_columns,M_columns>& y
		)
	{
		blitz::Range all = blitz::Range::all();
		blitz::TinyMatrix<T,N_rows,M_columns> z;
		for(int i=0; i < N_rows; ++i)
			for(int j=0; j < M_columns; ++j) {
				z(i,j) = 0;
				for(int k=0; k < N_columns; ++k)
					z(i,j) += x(i,k)*y(k,j);
			}
		
		return z;
	}

	template <class T_value, int N_elem> bool f_equal_to
		(const blitz::TinyVector<T_value,N_elem>& x1, 
		const blitz::TinyVector<T_value,N_elem>& x2)
	{
		return blitz_ext::is_each(x1,x2,std::equal_to<T_value>());
	}


	// return transposition of tiny matrix

	template <class T_num,int n_rows,int n_columns> 
	blitz::TinyMatrix<T_num,n_rows,n_columns>
		transpose(const blitz::TinyMatrix<T_num,n_rows,n_columns>& x) {

			blitz::TinyMatrix<T_num,n_rows,n_columns> y;
			for(int i_row = 0; i_row < n_rows; ++i_row)
				for(int i_col = 0; i_col < n_columns; ++i_col)
					y(i_col,i_row) = x(i_row,i_col);
			return y;

	}

	template <typename T_num1,typename T_num2,typename T_num3,int n_rows,int n_columns> 
	bool
		are_all_close(const blitz::TinyMatrix<T_num1,n_rows,n_columns>& x,
		const blitz::TinyMatrix<T_num2,n_rows,n_columns>& y,
		T_num3 rtol, T_num3 atol) {
			enum { n_elem = n_rows * n_columns };
			const T_num1 * p_x = x.data();
			const T_num2 * p_y = y.data();
			for(int i = 0; i < n_elem; ++i)
				if( ! PRODDL::Math::are_close(p_x[i],p_y[i],rtol,atol) )
					return false;
			return true;
	}

	template <typename T_num1,typename T_num2,int n_rows,int n_columns> 
	bool
		are_all_close(const blitz::TinyMatrix<T_num1,n_rows,n_columns>& x,
		const blitz::TinyMatrix<T_num2,n_rows,n_columns>& y) {
			return are_all_close(x,y,1e-8,1e-5);
	}

	template <typename T_num1,typename T_num2,typename T_num3,int n_elem> 
	bool
		are_all_close(const blitz::TinyVector<T_num1,n_elem>& x,
		const blitz::TinyVector<T_num2,n_elem>& y,
		T_num3 rtol, T_num3 atol) {
			const T_num1 * p_x = x.data();
			const T_num2 * p_y = y.data();
			for(int i = 0; i < n_elem; ++i)
				if( ! PRODDL::Math::are_close(p_x[i],p_y[i],rtol,atol) )
					return false;
			return true;
	}

	template <typename T_num1,typename T_num2,int n_elem> 
	bool
		are_all_close(const blitz::TinyVector<T_num1,n_elem>& x,
		const blitz::TinyVector<T_num2,n_elem>& y) {
			return are_all_close(x,y,1e-8,1e-5);
	}

} //namespace blitz_ext


#include <iostream>

namespace blitz {

	//   template<class T_num, int n_rows, int n_columns>
	//   std::ostream& operator<< (std::ostream& os, const TinyMatrix<T_num,n_rows,n_columns>& x)
	//   {
	//     os << n_rows << " x " << n_columns << std::endl << " [ ";
	//     for(int i_row = 0; i_row < n_rows; ++i_row) {
	//       os << " [ ";
	//       for(int i_col = 0; i_col < n_columns; ++i_col) {
	// 	os << std::setw(9) << x(i_row,i_col) << ' ';
	//       }
	//       os << " ] ";
	//     }
	//     os << " ] " << std::endl;
	//     return os;
	//   }

} // namespace blitz

#endif //AT_BZTINYVEC_H_
