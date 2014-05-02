//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_TYPES_H__
#define PRODDL_TYPES_H__

#include "PRODDL/Geom/traits.hpp"

#include "PRODDL/Geom/transformation.hpp"

#include "PRODDL/Grid/grid.hpp"

#include "PRODDL/Grid/project.hpp"

#include <complex>

namespace PRODDL { 

  template<typename T_num> struct Types : public Geom::SpaceTraits<T_num> {

    typedef Types<T_num> Self;
    
    typedef typename Self::Point3 Point;

    typedef typename Self::IntPoint3 IntPoint;

    typedef std::complex<T_num> Complex;

    // Grid of floating point numbers
    
    typedef Grid::Grid<3,T_num,T_num> Grid;

    typedef PRODDL::Grid::Grid<3,Complex,T_num> GridC;

    typedef PRODDL::Grid::Projector<T_num,T_num,3> Projector;

    typedef typename Self::VPoint3 Points;

    typedef typename common_types::num_vector_type<T_num>::Type Floats;

    typedef typename common_types::num_multiarray_type<T_num,2>::Type Matrix;

    typedef typename common_types::num_multiarray_type<int,2>::Type IntMatrix;

    typedef common_types::num_vector_type<int>::Type Ints;

    typedef common_types::num_vector_type<unsigned int>::Type Uints;

    typedef typename common_types::num_vector_type<Complex>::Type Complexes;

    typedef Geom::Translation<T_num> Translation;

    typedef Geom::Rotation<T_num> Rotation;

    typedef Geom::RotationTranslation<T_num> RotationTranslation;

    typedef typename common_types::num_vector_type<Translation>::Type Translations;

    typedef typename common_types::num_vector_type<Rotation>::Type Rotations;

    typedef typename common_types::num_vector_type<RotationTranslation>::Type RotationTranslations;

    struct TranValue {

      Translation tran;

      T_num value;

    };

    typedef typename common_types::num_vector_type<TranValue>::Type TranValues;

    struct RotTranValue {

      RotationTranslation tran;

      T_num value;

    };

    typedef typename common_types::num_vector_type<RotTranValue>::Type RotTranValues;

  }; // class Types<>

} //namespace PRODDL

#endif // PRODDL_TYPES_H__
