//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_BZ_CAST_H__
#define AT_BZ_CAST_H__

#include <blitz/array.h>

#include <exception>

#include <boost/static_assert.hpp>

namespace blitz {

  // Function that creates a view of blitz::Array<T_num1,n1> as blitz::Array<T_num2,n2>.
  // Primary purpose - view Array<TinyVector<T_num,n_tiny>,n_arr> as Array<T_num,n_arr+1>

  class bz_cast_error : public std::exception {
  protected:
    std::string m_msg;
  public:
    bz_cast_error(const std::string& msg) throw (): 
      m_msg(msg) {
    }
    virtual ~bz_cast_error() throw (){}
    virtual const char* what() const  throw () {
      return m_msg.c_str();
    }
  };

  namespace internal {
    template<typename T_num, int n_dim> void mergeTinyVecAndScalar(
								   const blitz::TinyVector<T_num,n_dim>& source_vec,
								   const T_num source_scal,
								   const int pos_scal,
								   blitz::TinyVector<T_num,n_dim+1>& dest_vec
								   )
    {
      for(int i=0; i < pos_scal; i++) dest_vec(i) = source_vec(i);
      dest_vec(pos_scal) = source_scal;
      for(int i=pos_scal+1; i < (n_dim+1); i++) dest_vec(i) = source_vec(i-1);
    }

    template<typename T_num, int n_dim> void removeFromTinyVec(
							       const blitz::TinyVector<T_num,n_dim>& source_vec,
							       const int pos,
							       blitz::TinyVector<T_num,n_dim-1>& dest_vec
							       )
    {
      for(int i=0; i < pos; i++) dest_vec(i) = source_vec(i);
      for(int i=pos+1; i < n_dim; i++) dest_vec(i-1) = source_vec(i);
    }

  template<int rank_o,int rank_n> struct dimension_folder {

    static void fold(const int stride_to_fold,
		     const TinyVector<int,rank_o>& stride_o, const TinyVector<int,rank_o>& shape_o, 
		     const GeneralArrayStorage<rank_o>& storage_o,
		     TinyVector<int,rank_n>& stride_n,  TinyVector<int,rank_n>& shape_n,
		     GeneralArrayStorage<rank_n>& storage_n) {

      BOOST_STATIC_ASSERT(rank_o >= rank_n);
      int stride_folded = 1;
      // new positions of old dimensions
      TinyVector<int,rank_o>& index;
      for(int k = 0; k < rank_o; k++) index[k] = k;
      int i = 0;
      // compute the new positions of old dimensions that will
      // be left and check that dimensions being deleted are
      // stored contiguously
      for( ; i < (rank_o - rank_n); i++) {
	int i_dim = storage_o.ordering()[i];
	if(stride_o[i_dim] != stride_folded * (storage_o.ascendingFlag()[i_dim] ? 1 : -1)) {
	  throw bz_cast_error("blitz::internal::folder::fold: Extents to be folded must be contiguous in memory.");
	}
	stride_folded *= shape_o[i_dim];
	for(int k = i_dim + 1; k < rank_o; k++) index[k] -= 1;
      }
      if( stride_folded != stride_to_fold ) {
	throw bz_cast_error("blitz::internal::folder::fold: Requested length to fold must be the product"\
			    " of (n_rank_old - n_rank_new) fastest extents.");
      }
      // copy the data of old dimensions that left
      for(int j=0; j < rank_n; j++, i++) {
	int ind_o = storage_o.ordering()[i];
	int ind_n = index[ind_o];
	stride_n[ind_n] = stride_o[ind_o]/stride_to_fold;
	shape_n [ind_n] = shape_o [ind_o];
	storage_n.base()[ind_n] = storage_o.base()[ind_o];
	storage_n.ascendingFlag()[ind_n] = storage_o.ascendingFlag()[ind_o];
	storage_n.ordering()[j] = ind_n;
      }
    }

    static void unfold(const int stride_to_unfold,
		       const TinyVector<int,rank_o>& stride_o, const TinyVector<int,rank_o>& shape_o, 
		     const GeneralArrayStorage<rank_o>& storage_o,
		     const TinyVector<int,rank_n-rank_o>& shape_d,
		     const GeneralArrayStorage<rank_n-rank_o>& storage_d,  
		     TinyVector<int,rank_n>& stride_n,  TinyVector<int,rank_n>& shape_n,
		     GeneralArrayStorage<rank_n>& storage_n) {


      BOOST_STATIC_ASSERT(rank_n > rank_o);
      // new positions of old dimensions
      TinyVector<int,rank_o>& index;
      for(int k = 0; k < rank_o; k++) index[k] = k;
      // copy data on additional dimensions to the new array descriptors and
      // find new positions for old dimensions
      int stride_unfolded = 1;
      int i=0;
      for( ; i < (rank_n - rank_o); i++) {
	int i_dim = storage_d.ordering()[i];
	// insertion of new dimension will move others one position up
	for(int k = i_dim; k < rank_o; k++) index[k] += 1;
	shape_n[i_dim] = shape_d[i];
	storage_n.base()[i_dim] = storage_d.base()[i];
	storage_n.ascendingFlag()[i_dim] = storage_d.ascendingFlag()[i];
	storage_n.ordering()[i] = i_dim;
	stride_n[i] = stride_unfolded * (storage_d.ascendingFlag()[i] ? 1 : -1);
	stride_unfolded *= shape_d[i];
      }
      if( stride_unfolded != stride_to_unfold ) {
	throw bz_cast_error("blitz::internal::folder::unfold: Shape of new dimensions must match the size of data type to unfold.");
      }
      // copy the data of old dimensions
      for(int j=0; j < rank_o; j++, i++) {
	int ind_o = storage_o.ordering()[i];
	int ind_n = index[ind_o];
	stride_n[ind_n] = stride_o[ind_o]*stride_to_unfold;
	shape_n [ind_n] = shape_o [ind_o];
	storage_n.base()[ind_n] = storage_o.base()[ind_o];
	storage_n.ascendingFlag()[ind_n] = storage_o.ascendingFlag()[ind_o];
	storage_n.ordering()[i] = ind_n;
      }
    }
  }; // struct folder


  } // namespace internal


  // Assuming that 'arr' is a multicomponent array in the same sence as in blitz::Array<>::extractComponent() method,
  // return new array of rank+1 that views old array's data as array of elements of components's type 'T_numtype_new'. 
  // New array inherits
  // storage format and bases from the old array. User must specify the position for new dimension as 'newDimPosition' parameter.
  // That means, indexing of new array will be done as A(ind_firstDimension,...,ind_newDimPosition,...). New dimension
  // will be considered stored ascending and with base zero.
  // New array is created from old array's data with memory policy 'neverDeleteData'. Thus, old array must exist during
  // lifetime of new array, or user must call new_array.reference(new_array.copy()) (looks like makeUnique() does not
  // work for arrays with 'nevereDeleteData') to create an independent copy.

  template<class T_num_o, class T_num_n, int rank_o> struct auto_rank {
    enum { sign_size = sizeof(T_num_n) > sizeof(T_num_o) ? 
	   1 : sizeof(T_num_n) < sizeof(T_num_o) ? -1 : 0 };
    enum { rank_n = rank_o - sign_size };
  };

  template<class T_num_o, class T_num_n, int rank_o, int rank_n, bool b_rank_d, bool b_size_d> class array_folder;

  template<class T_num_o, class T_num_n, int rank_o, int rank_n> 
  class array_folder<T_num_o,T_num_n,rank_o,rank_n,true,false> {
    typedef TinyVector<int,rank_n> Tiny_n;
    typedef TinyVector<int,rank_o> Tiny_o;
    typedef GeneralArrayStorage<rank_n> Storage_n;
    typedef GeneralArrayStorage<rank_o> Storage_o;
    typedef Array<T_num_n,rank_n> Array_n;
    typedef Array<T_num_o,rank_o> Array_o;

  };

  template<class T_numtype_new,class T_numtype, int N_rank>
  Array<T_numtype_new,N_rank+1> viewWithUnfoldedComponent(const Array<T_numtype,N_rank>& arr,
							  int newDimPosition = N_rank, int newDimBase = 0)
  {

    typedef TinyVector<diffType,N_rank+1> NewStride;
    typedef TinyVector<int,N_rank+1> NewShape;
    typedef GeneralArrayStorage<N_rank+1> NewStorage;

    typedef TinyVector<diffType,N_rank> OldStride;
    typedef TinyVector<int,N_rank> OldShape;
    typedef GeneralArrayStorage<N_rank> OldStorage;

    const std::size_t sizeT_numtype     = sizeof(T_numtype);
    const std::size_t sizeT_numtype_new = sizeof(T_numtype_new);

    // do compilation-time check that types are compatible
    BOOST_STATIC_ASSERT((sizeT_numtype >= sizeT_numtype_new) && (sizeT_numtype % sizeT_numtype_new == 0));

    const int numComponents = int(sizeT_numtype/sizeT_numtype_new);

    NewStride stride_new;
    internal::mergeTinyVecAndScalar(arr.stride(),diffType(0),newDimPosition,stride_new);
    stride_new *= numComponents;
    stride_new(newDimPosition) = 1;
    
    TinyVector<bool,N_rank> ascendingFlag_old;
    for(int i = 0; i < N_rank; i++) ascendingFlag_old(i) = arr.isRankStoredAscending(i);

    NewStorage storage_new;

    internal::mergeTinyVecAndScalar(ascendingFlag_old,true,newDimPosition,storage_new.ascendingFlag());

    for(int i = 0; i < N_rank; i++) {
      int dim_pos_old = arr.ordering(i);
      storage_new.ordering()(i+1) = (dim_pos_old < newDimPosition) ? dim_pos_old : (dim_pos_old + 1);
    }
    storage_new.ordering()(0) = newDimPosition;

    internal::mergeTinyVecAndScalar(arr.base(),newDimBase,newDimPosition,storage_new.base());

    const T_numtype_new* dataFirst_new = 
      reinterpret_cast<const T_numtype_new*>(arr.dataFirst());

    NewShape shape_new;
    internal::mergeTinyVecAndScalar(arr.shape(),numComponents,newDimPosition,shape_new);

    return Array<T_numtype_new,N_rank+1>(const_cast<T_numtype_new*>(dataFirst_new), 
					 shape_new, stride_new, 
					 neverDeleteData,storage_new);
  }


  template<class T_numtype_new,class T_numtype, int N_rank>
  Array<T_numtype_new,N_rank-1> viewWithFoldedComponent(const Array<T_numtype,N_rank>& arr)
  {

    typedef TinyVector<diffType,N_rank-1> NewStride;
    typedef TinyVector<int,N_rank-1> NewShape;
    typedef GeneralArrayStorage<N_rank-1> NewStorage;

    typedef TinyVector<diffType,N_rank> OldStride;
    typedef TinyVector<int,N_rank> OldShape;
    typedef GeneralArrayStorage<N_rank> OldStorage;

    const std::size_t sizeT_numtype     = sizeof(T_numtype);
    const std::size_t sizeT_numtype_new = sizeof(T_numtype_new);

    // do compilation-time test that type sizes are compatible
    BOOST_STATIC_ASSERT((sizeT_numtype <= sizeT_numtype_new) && (sizeT_numtype_new % sizeT_numtype == 0));

    const int numComponents = int(sizeT_numtype_new/sizeT_numtype);

    const int remDimPosition = arr.ordering(0);

    if( pow2(arr.stride(remDimPosition)) != 1 ) 
      throw bz_cast_error("blitz::viewWithFoldedComponent(): Array must have unit stride at the fastest changing dimension.");

    if( arr.extent(remDimPosition) != numComponents ) 
      throw bz_cast_error("blitz::viewWithFoldedComponent(): Extent of the fastest changing dimension must be equal to the number of components");

    NewStride stride_new;
    internal::removeFromTinyVec(arr.stride(),remDimPosition,stride_new);
    stride_new /= numComponents;
    
    TinyVector<bool,N_rank> ascendingFlag_old;
    for(int i = 0; i < N_rank; i++) ascendingFlag_old(i) = arr.isRankStoredAscending(i);

    NewStorage storage_new;

    internal::removeFromTinyVec(ascendingFlag_old,remDimPosition,storage_new.ascendingFlag());

    for(int i = 1; i < N_rank; i++) {
      int dim_pos_old = arr.ordering(i);
      storage_new.ordering()(i-1) = (dim_pos_old < remDimPosition) ? dim_pos_old : (dim_pos_old - 1);
    }

    internal::removeFromTinyVec(arr.base(),remDimPosition,storage_new.base());

    const T_numtype_new* dataFirst_new = 
      reinterpret_cast<const T_numtype_new*>(arr.dataFirst());

    NewShape shape_new;
    internal::removeFromTinyVec(arr.shape(),remDimPosition,shape_new);

    return Array<T_numtype_new,N_rank-1>(const_cast<T_numtype_new*>(dataFirst_new), 
					 shape_new, stride_new, 
					 neverDeleteData,storage_new);
  }

  template<class T_num, int N_rank_arr, int N_rank_tiny>
  Array<T_num,N_rank_arr+1> viewWithFlattenedComponent(const Array<TinyVector<T_num,N_rank_tiny>,N_rank_arr>& arr,
							  int newDimPosition = N_rank_arr, int newDimBase = 0)
  {
    return viewWithUnfoldedComponent<T_num>(arr,newDimPosition,newDimBase);
  }

} // namespace blitz

#endif // AT_BZ_CAST_H__
