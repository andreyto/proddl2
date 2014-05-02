/**
 * @file bob/core/array_utils.h
 * @date Tue Nov 8 15:34:31 2011 +0100
 * @author Andre Anjos <andre.anjos@idiap.ch>
 *
 * @brief Utilities for converting data to-from blitz::Arrays and other
 * goodies.
 *
 * Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BOB_CORE_ARRAY_UTILS_H 
#define BOB_CORE_ARRAY_UTILS_H

#include <blitz/array.h>
#include <stdint.h>
#include <stdexcept>
#include <boost/format.hpp>

#include <bob/core/array.h>
#include <bob/core/cast.h>

namespace bob { namespace core { namespace array {
  /**
   * @ingroup CORE_ARRAY
   * @{
   */

  /**
   * @brief Fills in shape and stride starting from a typeinfo object
   */
  template <int N> void set_shape_and_stride(const typeinfo& info,
      blitz::TinyVector<int,N>& shape, blitz::TinyVector<int,N>& stride) {
    for (int k=0; k<N; ++k) {
      shape[k] = info.shape[k];
      stride[k] = info.stride[k];
    }
  }


  /**
   * @brief Takes a data pointer and assumes it is a C-style array for the 
   * defined type. Creates a wrapper as a blitz::Array<T,N> with the same 
   * number of dimensions and type. Notice that the blitz::Array<> created 
   * will have its memory tied to the passed buffer. In other words you have 
   * to make sure that the buffer outlives the returned blitz::Array<>.
   */
  template <typename T, int N>
  blitz::Array<T,N> wrap(const interface& buf) {
    
    const typeinfo& type = buf.type();

    if (!buf.ptr()) throw std::runtime_error("empty buffer");

    if (type.dtype != bob::core::array::getElementType<T>()) {
      boost::format m("cannot efficiently retrieve blitz::Array<%s,%d> from buffer of type '%s'");
      m % stringize<T>() % N % type.str();
      throw std::runtime_error(m.str());
    }

    if (type.nd != N) {
      boost::format m("cannot retrieve blitz::Array<%s,%d> from buffer of type '%s'");
      m % stringize<T>() % N % type.str();
      throw std::runtime_error(m.str());
    }
          
    blitz::TinyVector<int,N> shape;
    blitz::TinyVector<int,N> stride;
    set_shape_and_stride(type, shape, stride);

    return blitz::Array<T,N>((T*)buf.ptr(), 
        shape, stride, blitz::neverDeleteData);
  }


  /**
   * @brief Takes a data pointer and assumes it is a C-style array for the 
   * defined type. Creates a copy as a blitz::Array<T,N> with the same number
   * of dimensions, but with a type as specified by you. If the type does not
   * match the type of the original C-style array, a cast will happen. 
   *
   * If a certain type cast is not supported. An appropriate exception will
   * be raised.
   */
  template <typename T, int N>
  blitz::Array<T,N> cast(const interface& buf) {

    const typeinfo& type = buf.type();

    if (type.nd != N) {
      boost::format m("cannot cast blitz::Array<%s,%d> from buffer of type '%s'");
      m % stringize<T>() % N % type.str();
      throw std::runtime_error(m.str());
    }

    switch (type.dtype) {
      case bob::core::array::t_bool: 
        return bob::core::array::cast<T>(wrap<bool,N>(buf));
      case bob::core::array::t_int8: 
        return bob::core::array::cast<T>(wrap<int8_t,N>(buf));
      case bob::core::array::t_int16: 
        return bob::core::array::cast<T>(wrap<int16_t,N>(buf));
      case bob::core::array::t_int32: 
        return bob::core::array::cast<T>(wrap<int32_t,N>(buf));
      case bob::core::array::t_int64: 
        return bob::core::array::cast<T>(wrap<int64_t,N>(buf));
      case bob::core::array::t_uint8: 
        return bob::core::array::cast<T>(wrap<uint8_t,N>(buf));
      case bob::core::array::t_uint16: 
        return bob::core::array::cast<T>(wrap<uint16_t,N>(buf));
      case bob::core::array::t_uint32: 
        return bob::core::array::cast<T>(wrap<uint32_t,N>(buf));
      case bob::core::array::t_uint64: 
        return bob::core::array::cast<T>(wrap<uint64_t,N>(buf));
      case bob::core::array::t_float32: 
        return bob::core::array::cast<T>(wrap<float,N>(buf));
      case bob::core::array::t_float64: 
        return bob::core::array::cast<T>(wrap<double,N>(buf));
      case bob::core::array::t_float128: 
        return bob::core::array::cast<T>(wrap<long double,N>(buf));
      case bob::core::array::t_complex64: 
        return bob::core::array::cast<T>(wrap<std::complex<float>,N>(buf));
      case bob::core::array::t_complex128: 
        return bob::core::array::cast<T>(wrap<std::complex<double>,N>(buf));
      case bob::core::array::t_complex256: 
        return bob::core::array::cast<T>(wrap<std::complex<long double>,N>(buf));
      default:
        break;
  }

    //if we get to this point, there is nothing much we can do...
    throw std::runtime_error("invalid type on blitz buffer array casting -- debug me");

  }

  /**
   * @}
   */
}}}

#endif /* BOB_CORE_ARRAY_UTILS_H */
