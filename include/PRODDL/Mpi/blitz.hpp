//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_MPI_BLITZ_H__
#define PRODDL_MPI_BLITZ_H__

// Methods and classes to pass Blitz datatype through MPI

#include "oompi.h"

#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include <blitz/array.h>

#include "PRODDL/exceptions.hpp"

#include "PRODDL/Blitz/bzarray.hpp"

#include "PRODDL/Mpi/tags.hpp"

#include "PRODDL/Mpi/mpi_logger.hpp"

namespace PRODDL {

  namespace Mpi {

    typedef char PrimitiveCastType;

    template<typename T_obj>
    OOMPI_Message
    messageContig(T_obj& obj,int tag) {
      return OOMPI_Message(reinterpret_cast<PrimitiveCastType*>(&obj),sizeof(T_obj),tag);
    }

    template<typename T_num, int n_dim>
    OOMPI_Message
    message(blitz::TinyVector<T_num,n_dim>& x,int tag = TAG_TINY_VECT) {
      return messageContig(x,tag);
    }

    template<typename T_num, int n_dim>
    OOMPI_Message
    message(blitz::Array<T_num,n_dim>& x,int tag = TAG_ARRAY_DATA) {
      if( ! (x.size() == 0 || x.isStorageContiguous()) ) {
#ifdef DBG_ENABLED
	x.dumpStructureInformation(dbg::out(dbg::error));
#endif
	dbg::assertion(dbg::error, DBG_ASSERTION(x.size() == 0 || x.isStorageContiguous()));
	throw not_supported_error("Mp::message(blitz::Array<T_num,n_dim>& x,int tag = TAG_ARRAY_DATA): Blitz array must be contiguous in memory.");
      }
      return OOMPI_Message(x.dataFirst(),x.size(),tag);
    }

    // Not working, reason unclear
//     template<typename T_num, int n_dim>
//     OOMPI_Message
//     messageArrayOfContigs(blitz::Array<T_num,n_dim>& x,int tag = TAG_ARRAY_DATA) {
//       if( ! x.isStorageContiguous() )
// 	throw not_supported_error("Mpi: Blitz array must be contiguous in memory.");
//       OOMPI_Message msgElement(message(*x.dataFirst(),tag));
//       return OOMPI_Message(OOMPI_Datatype(msgElement),x.dataFirst(),x.size(),tag);
//     }


    template<typename T_num, int n_dim>
    OOMPI_Message
    messageArrayOfContigs(blitz::Array<T_num,n_dim>& x,int tag = TAG_ARRAY_DATA) {
      if( ! (x.size() == 0 || x.isStorageContiguous()) ) {
#ifdef DBG_ENABLED
	x.dumpStructureInformation(dbg::out(dbg::error));
#endif
	dbg::assertion(dbg::error, DBG_ASSERTION(x.size() == 0 || x.isStorageContiguous()));
	throw not_supported_error("Mpi::messageArrayOfContigs(blitz::Array<T_num,n_dim>& x,int tag = TAG_ARRAY_DATA) : Blitz array must be contiguous in memory.");
      }
      return OOMPI_Message(reinterpret_cast<PrimitiveCastType*>(x.dataFirst()),x.size()*sizeof(T_num),tag);
    }


    template<typename T_num, int n_dim_arr, int n_dim_tiny>
    OOMPI_Message
    message(blitz::Array<blitz::TinyVector<T_num,n_dim_tiny>, n_dim_arr>& x,int tag = TAG_ARRAY_DATA) {
      return messageArrayOfContigs(x,tag);      
    }

    template<typename T_val,int n_dim>
    int elements(OOMPI_Status& status, OOMPI_Message& message, 
		 const blitz::Array<T_val,n_dim>& dummy) {

      return status.Get_elements(OOMPI_Datatype(message.Get_type()))/(sizeof(T_val)/sizeof(PrimitiveCastType));

    }


    template<typename T_num, int n_dim>
    class ArrayMessenger {

    public:

      typedef blitz::Array<T_num,n_dim> ArrayType;
      typedef ArrayStruct<n_dim> ArrayStructType;

    public:
      
      ArrayMessenger(ArrayType& _arr):
	p_arr(&_arr),
	arr_struct(_arr)
      {}


      OOMPI_Message
      messageStruct() { 
	return messageContig(arr_struct,TAG_ARRAY_STRUCT);
      }

      OOMPI_Message
      messageData(bool do_validate = true) {

	if( do_validate ) {
	  if( structChanged() )
	    rebuildArray();
	}

	return message(*p_arr);

      }

    protected:


      bool structChanged() {

	// if anything in structure received in messages is different from 
	// a structure of array referenced by this object

	return ! arr_struct.isSameStruct(*p_arr);

      }

      void rebuildArray() {
	
	p_arr->reference(arr_struct.createArray(T_num()));

      }

    protected:
      
      ArrayType * p_arr;

      ArrayStructType arr_struct;

    }; // class ArrayMessenger


  } } // namespace PRODDL::Mpi

#endif // PRODDL_MPI_BLITZ_H__
