//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_GEOM_MOVE_H__
#define AT_GEOM_MOVE_H__

#include "PRODDL/Geom/traits.hpp"
#include "PRODDL/Geom/geom_except.hpp"
#include "PRODDL/Common/math.hpp"

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Geom/euler.hpp"

#include "PRODDL/Geom/transformation.hpp"

#include "PRODDL/External/nr/nr_arr.hpp"
#include "PRODDL/External/nr/nr.hpp"

#include <algorithm> // min & max

namespace PRODDL { namespace Geom { namespace Transforms {


  //Perform coordinate transformation using rotation matrix mRot, translation vector vTran and geom. center
  //vector vCM_in on vectors in array avIn and put result in array avOut.
  //avOut is allowed to point to the same data as avIn.

  template<typename T_num>
  inline 
  void applyRotTran(const typename SpaceTraits<T_num>::Matrix3x3& mRot,
		    const typename SpaceTraits<T_num>::Point3& vTran,
		    const typename SpaceTraits<T_num>::Point3& vCM_in,
		    const typename SpaceTraits<T_num>::VPoint3& avIn,
		    typename SpaceTraits<T_num>::VPoint3& avOut)
  {
    typename SpaceTraits<T_num>::Point3 vTr = vCM_in - blitz::product(mRot,vCM_in) + vTran;
    if(avIn.size() != avOut.size()) {
      throw transforms_size_error("applyRotTran(): avIn and avOut array sizes do not match");
    }
    int n = avIn.rows();
    for(int i = 0; i < n; i++)
      {
	avOut(i) = blitz::product(mRot,avIn(i)) + vTr;
      }
  }


  template<typename T_num>
  class RigidBody {

  public:

    typedef typename Geom::SpaceTraits<T_num>::Point3 Point;
    typedef typename Geom::SpaceTraits<T_num>::Matrix3x3 Matrix3;
    typedef typename Geom::SpaceTraits<T_num>::VPoint3 Points;
    typedef typename Geom::SpaceTraits<T_num>::Point3Pair PointPair;

    typedef typename common_types::num_vector_type<T_num>::Type fvect;
    typedef typename common_types::point_type<bool,3>::Type BoolPoint;

    typedef typename Geom::TransformationTraits<T_num>::VRotationTranslation VRotationTranslation;
    typedef Geom::Rotation<T_num> Rotation;
    typedef Geom::Translation<T_num> Translation;
    typedef Geom::RotationTranslation<T_num> RotationTranslation;


    // Convert Cartesian gradients to rigid body gradients.
    // Take positions 'xyzCoords' and gradients 'xyzGrads' in Cartesian coords for a set of points, 
    // as well as position 'rbCoords' in rigid body coords (entire set is considered to be 
    // a rigid body). Set output parameter 'rbGrads' to gradients in rigid body coords 
    // (first element of rbGrads is translational gradient, second - rotational).
    // Requirements: 'rbCoords(1)' angles are assumed to define rotation around the point 'rbCoords(0)'.

    static
    void xyzGradToRigid(const Points& xyzCoords,
			const Points& xyzGrads,
			const PointPair& rbCoords,
			PointPair& rbGrads)
    {


      int n_points = xyzCoords.size();

      rbGrads(0) = 0;
      rbGrads(1) = 0;

      //    compute the rigid body gradient components for each group
      //
      T_num phi = rbCoords(1)(0);
      T_num theta = rbCoords(1)(2); // Tinker's is swapped with psi compared to our convention
      T_num cphi = std::cos(phi);
      T_num sphi = std::sin(phi);
      T_num ctheta = std::cos(theta);
      T_num stheta = std::sin(theta);
      //
      //    get unit vectors along the phi, theta and psi rotation axes
      //
      Point ephi, etheta, epsi;
      ephi(0) = 0.0;
      ephi(1) = 0.0;
      ephi(2) = 1.0;
      etheta(0) = -sphi;
      etheta(1) = cphi;
      etheta(2) = 0.0;
      epsi(0) = ctheta * cphi;
      epsi(1) = ctheta * sphi;
      epsi(2) = -stheta;
      //
      //    first, get the rigid body gradients for translations
      //
      for(int i=0; i < n_points; i++) {
	rbGrads(0) += xyzGrads(i);
      }
      //
      //    accumulate the moment arm along each axis of rotation
      //
      Point tau;
      tau = 0;
      for(int i=0; i < n_points; i++) {
	Point xyzTerm = xyzCoords(i) - rbCoords(0);
	const Point& xyzGrad = xyzGrads(i);
	tau(0) += xyzTerm(1)*xyzGrad(2) - xyzTerm(2)*xyzGrad(1);
	tau(1) += xyzTerm(2)*xyzGrad(0) - xyzTerm(0)*xyzGrad(2);
	tau(2) += xyzTerm(0)*xyzGrad(1) - xyzTerm(1)*xyzGrad(0);
      }
      //
      //    now, set the rigid body gradients for rotations
      //
      for(int j=0; j < 3; j++) {
	rbGrads(1) += Point(ephi(j),etheta(j),epsi(j))*tau(j);
      }
      //
      // swap the last two angle components because of a different storage convention
      // compared to Tinker
      //
      T_num dummy = rbGrads(1)(2);
      rbGrads(1)(2) = rbGrads(1)(1);
      rbGrads(1)(1) = dummy;

#if ATDEBUG_LEVEL > 8

      ATOUTVAR(rbGrads); ATOUTENDL();

#endif

    }


    /*
      c
      c
      c     ##############################################################
      c     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ##
      c     ##                   All Rights Reserved                    ##
      c     ##############################################################
      c
      c     ###############################################################
      c     ##                                                           ##
      c     ##  subroutine standardOrientation  --  rigid body reference coordinates  ##
      c     ##                                                           ##
      c     ###############################################################
      c
      c
      c     "standardOrientation" computes a set of reference Cartesian coordinates
      c     in standard orientation for each rigid body atom group
      c
      c
    */

    static
    void standardOrientation(const Points& xyz, const fvect& mass, Points& xyzrb, PointPair& rbc) {

      /*
	c
	c
	c     use current coordinates as default reference coordinates
	c
      */

      //xyzrb = xyz;

      /*
	c
	c     compute the rigid body coordinates for each atom group
	c
      */
      xyzrigid(xyz,mass,rbc);

      /*
	c
	c     translate and rotate each atom group into inertial frame
	c
      */

      //(RotationTranslation(rbc).inverse())(xyzrb);

      /*
	c
	c     get the center of mass and Euler angles for each group     
	c
      */

	Point xyzcm = rbc(0);

	T_num       phi = rbc(1)(0);
	T_num       theta = rbc(1)(2);
	T_num       psi = rbc(1)(1);
	T_num       cphi = std::cos(phi);
	T_num       sphi = std::sin(phi);
	T_num       ctheta = std::cos(theta);
	T_num       stheta = std::sin(theta);
	T_num       cpsi = std::cos(psi);
	T_num       spsi = std::sin(psi);

      /*
	c
	c     construct the rotation matrix from Euler angle values
	c
      */

	Matrix3 a;
	a(0,0) = ctheta * cphi;
	a(1,0) = spsi*stheta*cphi - cpsi*sphi;
	a(2,0) = cpsi*stheta*cphi + spsi*sphi;
	a(0,1) = ctheta * sphi;
	a(1,1) = spsi*stheta*sphi + cpsi*cphi;
	a(2,1) = cpsi*stheta*sphi - spsi*cphi;
	a(0,2) = -stheta;
	a(1,2) = ctheta * spsi;
	a(2,2) = ctheta * cpsi;

      /*
	c
	c     translate and rotate each atom group into inertial frame
	c
      */

	for(int k=0;k<xyz.size();k++) {
	 T_num xterm = xyz(k)(0) - xyzcm(0);
	 T_num yterm = xyz(k)(1) - xyzcm(1);
	 T_num zterm = xyz(k)(2) - xyzcm(2);
	 xyzrb(k)(0) = a(0,0)*xterm + a(0,1)*yterm + a(0,2)*zterm;
	 xyzrb(k)(1) = a(1,0)*xterm + a(1,1)*yterm + a(1,2)*zterm;
	 xyzrb(k)(2) = a(2,0)*xterm + a(2,1)*yterm + a(2,2)*zterm;
	}

    }

    /*
      c
      c
      c     #################################################################
      c     ##                                                             ##
      c     ##  subroutine xyzrigid  --  determine rigid body coordinates  ##
      c     ##                                                             ##
      c     #################################################################
      c
      c
      c     "xyzrigid" computes the center of mass and Euler angle rigid
      c     body coordinates for each atom group in the system
      c
      c     literature reference:
      c
      c     Herbert Goldstein, "Classical Mechanics, 2nd Edition",
      c     Addison-Wesley, Reading, MA, 1980; see the Euler angle
      c     xyz convention in Appendix B
      c
      c
    */

    static
    void xyzrigid(const Points& xyz, const fvect& mass, PointPair& rbc)
    {

      /*
	c
	c     compute the position of the group center of mass
	c
      */

      Point xyzcm = blitz::mean(xyz*mass);

      /*
	c
	c     compute and then diagonalize the inertia tensor
	c
      */

      T_num    xx = 0.0;
      T_num    xy = 0.0;
      T_num    xz = 0.0;
      T_num    yy = 0.0;
      T_num    yz = 0.0;
      T_num    zz = 0.0;

      for(int k = 0; k < xyz.size(); k++) {

	T_num weigh = mass(k);
	T_num xterm = xyz(k)(0) - xyzcm(0);
	T_num yterm = xyz(k)(1) - xyzcm(1);
	T_num zterm = xyz(k)(2) - xyzcm(2);
	xx = xx + xterm*xterm*weigh;
	xy = xy + xterm*yterm*weigh;
	xz = xz + xterm*zterm*weigh;
	yy = yy + yterm*yterm*weigh;
	yz = yz + yterm*zterm*weigh;
	zz = zz + zterm*zterm*weigh;

      }

      Matrix3 tensor;
      Matrix3 vec;

      tensor(0,0) = yy + zz;
      tensor(1,0) = -xy;
      tensor(2,0) = -xz;
      tensor(0,1) = -xy;
      tensor(1,1) = xx + zz;
      tensor(2,1) = -yz;
      tensor(0,2) = -xz;
      tensor(1,2) = -yz;
      tensor(2,2) = xx + yy;

      //         call jacobi (3,3,tensor,moment,vec,work1,work2)
      Point moment;
      int nrot;
      PRODDL::nr::MatrixAdaptor<T_num> tensor_nrc(tensor.dataFirst(),1,3,1,3), vec_nrc(vec.dataFirst(),1,3,1,3);
      //jacoby diagonalization
      PRODDL::nr::jacobi(tensor_nrc.rowPointers(),3,moment.dataFirst()-1,vec_nrc.rowPointers(),&nrot); 
      //sort eigenvalues in ascending order, like
      //Tinker's jacobi does
      PRODDL::nr::eigsrt(moment.dataFirst()-1,vec_nrc.rowPointers(),3,true); 

      /*
	c
	c     select the direction for each principle moment axis
	c
      */

      for(int m = 0; m < 2; m++) {
	for(int k = 0; k < xyz.size(); k++) {
	
	  T_num xterm = vec(0,m) * (xyz(k)(0)-xyzcm(0));
	  T_num yterm = vec(1,m) * (xyz(k)(1)-xyzcm(1));
	  T_num zterm = vec(2,m) * (xyz(k)(2)-xyzcm(2));
	  T_num dot = xterm + yterm + zterm;
	  if (dot < 0.0) {
	    vec(0,m) = -vec(0,m);
	    vec(1,m) = -vec(1,m);
	    vec(2,m) = -vec(2,m);
	  }
	  if (dot != 0.0)  break;
	}
      }

      /*
	c
	c     moment axes must give a right-handed coordinate system
	c
      */
      T_num xterm = vec(0,0) * (vec(1,1)*vec(2,2)-vec(1,2)*vec(2,1));
      T_num yterm = vec(1,0) * (vec(0,2)*vec(2,1)-vec(0,1)*vec(2,2));
      T_num zterm = vec(2,1) * (vec(0,1)*vec(1,2)-vec(0,2)*vec(1,1));
      T_num dot = xterm + yterm + zterm;
      if (dot < 0.0) {
	for(int j = 0; j < 3; j++) {
	  vec(j,2) = -vec(j,2);
	}
      }

#if ATDEBUG_LEVEL > 8
      ATOUTVAR(vec); ATOUTENDL();
      {
	Point 
	  a1(vec(0,0),vec(1,0),vec(2,0)),
	  a2(vec(0,1),vec(1,1),vec(2,1)),
	  a3(vec(0,2),vec(1,2),vec(2,2));
	ATOUTVAR(blitz::dot(a1,a2)); ATOUTVAR(blitz::dot(a1,a3)); ATOUTVAR(blitz::dot(a2,a3)); ATOUTENDL();
      }
      {
	Point 
	  a1(vec(0,0),vec(0,1),vec(0,2)),
	  a2(vec(1,0),vec(1,1),vec(1,2)),
	  a3(vec(2,0),vec(2,1),vec(2,2));
	ATOUTVAR(blitz::dot(a1,a2)); ATOUTVAR(blitz::dot(a1,a3)); ATOUTVAR(blitz::dot(a2,a3)); ATOUTENDL();
      }
#endif

      /*
	c
	c     principal moment axes form rows of Euler rotation matrix
	c
      */
    
      Matrix3 a;
      for(int k = 0; k < 3; k++)
	for(int j = 0; j < 3; j++)
	  a(k,j) = vec(j,k);

      /*
	c
	c     compute Euler angles consistent with the rotation matrix
	c
      */

      roteuler(a,rbc(1));
      //rbc(1) = Eu_xyz::rotMatrToEulerAng<T_num>(a);

      /*
	c
	c     set the rigid body coordinates for each atom group
	c
      */

      rbc(0) = xyzcm;

    }

    /*
      c
      c
      c     #################################################################
      c     ##                                                             ##
      c     ##  subroutine roteuler  --  rotation matrix to Euler angles   ##
      c     ##                                                             ##
      c     #################################################################
      c
      c
      c     "roteuler" computes a set of Euler angle values consistent
      c     with an input rotation matrix
      c
      c
    */

    static
    void roteuler (const Matrix3& a, Point& ang) {
      T_num &phi = ang(0), &psi = ang(1), &theta = ang(2);

      /*
	c
	c
	c     set the tolerance for Euler angles and rotation elements
	c
      */

      T_num eps = 1.0e-8;

      /*
	c
	c     get a trial value of theta from a single rotation element
	c
      */

      theta = std::asin(std::min(1.0,std::max(-1.0,-a(0,2))));
      T_num ctheta = std::cos(theta);
      T_num stheta = -a(0,2);

      /*
	c
	c     set the phi/psi difference when theta is either 90 or -90
	c
      */

      if (std::abs(ctheta) <= eps) {
	phi = 0.0;
	if (std::abs(a(2,0)) < eps) {
	  psi = std::asin(std::min(1.0,std::max(-1.0,-a(1,0)/a(0,2))));
	}
	else if (std::abs(a(1,0)) < eps) {
	  psi = std::acos(std::min(1.0,std::max(-1.0,-a(2,0)/a(0,2))));
	}
	else {
	  psi = std::atan(a(1,0)/a(2,0));
	}
      }

      /*
	c
	c     set the phi and psi values for all other theta values
	c
      */

      else {
	if (std::abs(a(0,0)) < eps) {
	  phi = std::asin(std::min(1.0,std::max(-1.0,a(0,1)/ctheta)));
	}
	else if (std::abs(a(0,1)) < eps) {
	  phi = std::acos(std::min(1.0,std::max(-1.0,a(0,0)/ctheta)));
	}
	else {
	  phi = std::atan(a(0,1)/a(0,0));
	}
	if (std::abs(a(2,2)) < eps) {
	  psi = std::asin(std::min(1.0,std::max(-1.0,a(1,2)/ctheta)));
	}
	else if (std::abs(a(1,2)) < eps) {
	  psi = std::acos(std::min(1.0,std::max(-1.0,a(2,2)/ctheta)));
	}
	else {
	  psi = std::atan(a(1,2)/a(2,2));
	}
      }

      /*
	c
	c     find sine and cosine of the trial phi and psi values
	c
      */

      T_num cphi = std::cos(phi);
      T_num sphi = std::sin(phi);
      T_num cpsi = std::cos(psi);
      T_num spsi = std::sin(psi);

      /*
	c
	c     reconstruct the diagonal of the rotation matrix
	c
      */

      Point b;
      b(0) = ctheta * cphi;
      b(1) = spsi*stheta*sphi + cpsi*cphi;
      b(2) = ctheta * cpsi;

      /*
	c
	c     compare the correct matrix diagonal to rebuilt diagonal
	c
      */


      BoolPoint err;
      for(int i = 0; i < 3; i++) {
	err(i) = false;
	if (std::abs(a(i,i)-b(i)) > eps)  err(i) = true;
      }

      /*
	c
	c     alter Euler angles to get correct rotation matrix values
	c
      */

      T_num pi = T_num(4.0)*std::atan(T_num(1.0));

      if (err(0) && err(1))  phi = phi - Math::sign(pi,phi);
      if (err(0) && err(2))  theta = -theta + Math::sign(pi,theta);
      if (err(1) && err(2))  psi = psi - Math::sign(pi,psi);

    }


    /*
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine rigidxyz  --  rigid body to Cartesian coords  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "rigidxyz" computes Cartesian coordinates for a rigid body
c     group via rotation and translation of reference coordinates  
c
c     literature reference:
c
c     Herbert Goldstein, "Classical Mechanics, 2nd Edition",
c     Addison-Wesley, Reading, MA, 1980; see the Euler angle
c     xyz convention in Appendix B
c
c
    */

    static
    void rigidxyz(const Points& xyzrb, const PointPair& rbc, Points& xyz)
    {

      /*
c
c
c     get the center of mass and Euler angles for each group     
c
      */

	Point xyzcm = rbc(0);

	T_num       phi = rbc(1)(0);
	T_num       theta = rbc(1)(2);
	T_num       psi = rbc(1)(1);


	T_num       cphi = std::cos(phi);
	T_num       sphi = std::sin(phi);
	T_num       ctheta = std::cos(theta);
	T_num       stheta = std::sin(theta);
	T_num       cpsi = std::cos(psi);
	T_num       spsi = std::sin(psi);

	/*
c
c     construct the rotation matrix from Euler angle values
c
	*/

	Matrix3 a;
	a(0,0) = ctheta * cphi;
	a(1,0) = spsi*stheta*cphi - cpsi*sphi;
	a(2,0) = cpsi*stheta*cphi + spsi*sphi;
	a(0,1) = ctheta * sphi;
	a(1,1) = spsi*stheta*sphi + cpsi*cphi;
	a(2,1) = cpsi*stheta*sphi - spsi*cphi;
	a(0,2) = -stheta;
	a(1,2) = ctheta * spsi;
	a(2,2) = ctheta * cpsi;

	/*
c
c     rotate and translate reference coordinates into global frame
c
	*/

	for(int k=0;k<xyz.size();k++) {
	  T_num xterm = xyzrb(k)(0);
	  T_num yterm = xyzrb(k)(1);
	  T_num zterm = xyzrb(k)(2);
	  xyz(k)(0) = a(0,0)*xterm + a(1,0)*yterm + a(2,0)*zterm + xyzcm(0);
	  xyz(k)(1) = a(0,1)*xterm + a(1,1)*yterm + a(2,1)*zterm + xyzcm(1);
	  xyz(k)(2) = a(0,2)*xterm + a(1,2)*yterm + a(2,2)*zterm + xyzcm(2);
	}
    }

    static
    RotationTranslation transformation(const PointPair& rbc) {

	Point xyzcm = rbc(0);

	T_num       phi = rbc(1)(0);
	T_num       theta = rbc(1)(2);
	T_num       psi = rbc(1)(1);


	T_num       cphi = std::cos(phi);
	T_num       sphi = std::sin(phi);
	T_num       ctheta = std::cos(theta);
	T_num       stheta = std::sin(theta);
	T_num       cpsi = std::cos(psi);
	T_num       spsi = std::sin(psi);

	/*
c
c     construct the rotation matrix from Euler angle values
c
	*/

	Matrix3 a;
	a(0,0) = ctheta * cphi;
	a(1,0) = spsi*stheta*cphi - cpsi*sphi;
	a(2,0) = cpsi*stheta*cphi + spsi*sphi;
	a(0,1) = ctheta * sphi;
	a(1,1) = spsi*stheta*sphi + cpsi*cphi;
	a(2,1) = cpsi*stheta*sphi - spsi*cphi;
	a(0,2) = -stheta;
	a(1,2) = ctheta * spsi;
	a(2,2) = ctheta * cpsi;
	
	return RotationTranslation(Math::transpose(a),xyzcm);
    }
    

  }; // class RigidBody

} // namespace Transforms


} } //namespace PRODDL::Geom

#endif // AT_GEOM_MOVE_H__
