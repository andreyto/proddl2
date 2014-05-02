//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_SYMMETRY_H__
#define PRODDL_SYMMETRY_H__

#include "PRODDL/Geom/transformation.hpp"

namespace PRODDL { namespace Geom { 


  template<typename T_num>
  class SymmetryCheckerCn {

  public:

    typedef typename SpaceTraits<T_num>::Point3 Point;
    typedef RotationTranslation<T_num> RotationTranslationT;

    SymmetryCheckerCn(int nSymm, T_num radius) :
      m_nSymm(nSymm),
      unit1(radius,0,0), 
      unit2(0,radius,0), 
      unit3(0,0,radius), 
      unit4(0,0,0) 
    {}


    T_num operator() (RotationTranslationT transform) const {

      RotationTranslationT transformN = transform;

      for(int i_symm = 1; i_symm < m_nSymm; i_symm++) {
	transformN = transform * transformN;
      }

      T_num rmsd = 
	Math::dotSelf(transformN(unit1) - unit1) +
	Math::dotSelf(transformN(unit2) - unit2) +
	Math::dotSelf(transformN(unit3) - unit3) +
	Math::dotSelf(transformN(unit4) - unit4);

      return std::sqrt(rmsd/4);

    }

    int getNSymm() const {
      return m_nSymm;
    }

  protected:

    int m_nSymm;

    Point unit1, unit2, unit3, unit4;

  }; // class SymmetryCheckerCn


} } // namespace PRODDL { namespace Geom {

#endif // PRODDL_SYMMETRY_H__
