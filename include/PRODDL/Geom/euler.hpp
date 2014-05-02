//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_GEOM_EULER_H__
#define AT_GEOM_EULER_H__

#include "PRODDL/Geom/traits.hpp"
#include "PRODDL/Geom/geom_except.hpp"
#include "PRODDL/Common/math.hpp"


namespace PRODDL { namespace Geom { namespace Transforms {


  //Separate namespace is used for each euler angles convention:

  //Invariants:
  //eulerAngToRotMatr(rotMatrToEulerAng(eulerAngToRotMatr(angles))) == eulerAngToRotMatr(angles)

  //Warning: when creating and using these functions, attention should be paid to the
  //Euler angles convention (XYZ, X, ZYX, fixed frame, body frame etc).
  //In particular, code in MagicSoftware defines XYZ and others differently from GRAMM
  //and Mathematica. Tinker stores in arrays (Phi,Theta,Psi) in this order, otherwise
  //convention is the same as here. In other words, to adopt Tinker algorithms,
  //the second and third angle components must be swapped.

  //Eu_xyz - yown[-Pi,Pi], roll[-Pi,Pi], pitch[-Pi/2,Pi/2] or Phi, Psi, Theta - the ones used in gramm

  namespace Eu_xyz
  {

    //create rotation matrix from Euler angles (angles must be expressed in radians):
    template<typename T_num> 
    typename SpaceTraits<T_num>::Matrix3x3 
    eulerAngToRotMatr(const typename SpaceTraits<T_num>::Point3& angles)
    {
      T_num a = angles[0], b = angles[1], g = angles[2];
      T_num sa = sin(a), ca = cos(a),
	sb = sin(b), cb = cos(b),
	sg = sin(g), cg = cos(g);
      typename SpaceTraits<T_num>::Matrix3x3 m;
      m =
	cg*ca,          cg*sa,          -sg,
	sb*sg*ca-cb*sa, sb*sg*sa+cb*ca, cg*sb,
	cb*sg*ca+sb*sa, cb*sg*sa-sb*ca, cg*cb;
      return m;
    }

    //Create XYZ Euler angles from rotation matrix (angles will be expressed in radians)
    //Arguments: Matrix3x3 - rotation matrix [Input], bool& isResultUnique - false if
    //returned angles are not unique [Output]
    //Return: Point3 - Euler angles
    template<typename T_num> 
    typename SpaceTraits<T_num>::Point3
    rotMatrToEulerAng (const typename SpaceTraits<T_num>::Matrix3x3& m, bool& isResultUnique)
    {
      // rot =  cy*cz          -cy*sz           sy
      //        cz*sx*sy+cx*sz  cx*cz-sx*sy*sz -cy*sx
      //       -cx*cz*sy+sx*sz  cz*sx+cx*sy*sz  cx*cy

      T_num HALF_PI = T_num(4.0)*std::atan(T_num(1.0))/T_num(2.);

      typename SpaceTraits<T_num>::Point3 angles;
      if ( m(0,2) < 1.0 )
	{
	  if ( m(0,2) > -1.0 )
	    {
	      angles(0) = std::atan2(-m(1,2),m(2,2));
	      angles(1) = std::asin(m(0,2));
	      angles(2) = std::atan2(-m(0,1),m(0,0));
	      isResultUnique =  true;
	    }
	  else
	    {
	      // WARNING.  Not unique.  XA - ZA = atan2(r10,r11)
	      angles(0) = -std::atan2(m(1,0),m(1,1));
	      angles(1) = -HALF_PI;
	      angles(2) = 0.0;
	      isResultUnique = false;
	    }
	}
      else
	{
	  // WARNING.  Not unique.  XA + ZA = atan2(r10,r11)
	  angles(0) = std::atan2(m(1,0),m(1,1));
	  angles(1) = HALF_PI;
	  angles(2) = 0.0;
	  isResultUnique = false;
	}
    
      //Swap angles to convert from Magic's -Y-P-R to our PRY
      typename SpaceTraits<T_num>::Point3 our_angles;
      our_angles(0) = - angles(2); 
      our_angles(1) = - angles(0); 
      our_angles(2) = - angles(1);
      return our_angles;
    }

    template<typename T_num> 
    typename SpaceTraits<T_num>::Point3
    rotMatrToEulerAng (const typename SpaceTraits<T_num>::Matrix3x3& m) {
      bool dummy;
      return rotMatrToEulerAng<T_num>(m,dummy);
    }

  } //end namespace Eu_xyz


}}} // namespace


#endif // PRODDL_GEOM_EULER_H__
