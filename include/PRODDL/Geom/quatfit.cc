
/*
   c     "quatfit" uses a quaternion based method to achieve the best
   c     fit superposition of two sets of coordinates
   c
   c     literature reference:
   c
   c     S. J. Kearsley, "An Algorithm for the Simultaneous Superposition
   c     of a Structural Series", Journal of Computational Chemistry,
   c     11, 1187-1192 (1990)
   c
   c     adapted from Tinker quatfit.f (J.W.Ponder) (1992)
   c     adapted from an original program written by David J. Heisterberg,
   c     Ohio Supercomputer Center, Columbus, OH
*/
#include "DockTK/Geom/quatfit.hpp"

#include "DockTK/External/nr/nr_arr.hpp"
#include "DockTK/External/nr/nr.hpp"

#include "DockTK/Common/math.hpp"

namespace DockTK { namespace Geom { namespace Align {

template<typename T_num,class PointsInputIterator1,class PointsInputIterator2,class WeightsInputIterator>
typename QuatFitTraits<T_num>::Matr3 
quatFit(PointsInputIterator1 points1_first,
	PointsInputIterator1 points1_last,
	PointsInputIterator2 points2_first,
	WeightsInputIterator weights_first) {

  typedef T_num mfl;
  typedef blitz::TinyVector<T_num,3> vect;
  typedef blitz::TinyMatrix<T_num,4,4> matrf4;
  typedef blitz::TinyVector<T_num,4> vect4;
   mfl xxyx = 0;
   mfl xxyy = 0;
   mfl xxyz = 0;
   mfl xyyx = 0;
   mfl xyyy = 0;
   mfl xyyz = 0;
   mfl xzyx = 0;
   mfl xzyy = 0;
   mfl xzyz = 0;
   
   mfl w_sum = 0;
   for( ; points1_first != points1_last; ++points1_first, ++points2_first, ++weights_first)
   {
      mfl w = *weights_first;
      vect v1 = *points1_first;
      vect v2 = *points2_first;
      xxyx +=w*v1(0)*v2(0);
      xxyy +=w*v1(1)*v2(0);
      xxyz +=w*v1(2)*v2(0);
      xyyx +=w*v1(0)*v2(1);
      xyyy +=w*v1(1)*v2(1);
      xyyz +=w*v1(2)*v2(1);
      xzyx +=w*v1(0)*v2(2);
      xzyy +=w*v1(1)*v2(2);
      xzyz +=w*v1(2)*v2(2);

      w_sum += w;
   }

   if( w_sum != 0 ) {
     xxyx /= w_sum;
     xxyy /= w_sum;
     xxyz /= w_sum;
     xyyx /= w_sum;
     xyyy /= w_sum;
     xyyz /= w_sum;
     xzyx /= w_sum;
     xzyy /= w_sum;
     xzyz /= w_sum;
   }

   matrf4 c,v;
   c(0,0) = xxyx + xyyy + xzyz;
   c(0,1) = xzyy - xyyz;
   c(1,1) = xxyx - xyyy - xzyz;
   c(0,2) = xxyz - xzyx;
   c(1,2) = xxyy + xyyx;
   c(2,2) = xyyy - xzyz - xxyx;
   c(0,3) = xyyx - xxyy;
   c(1,3) = xzyx + xxyz;
   c(2,3) = xyyz + xzyy;
   c(3,3) = xzyz - xxyx - xyyy;

   int nrot;
   vect4 d;
   DockTK::nr::MatrixAdaptor<T_num> c_nrc(c.dataFirst(),1,4,1,4), v_nrc(v.dataFirst(),1,4,1,4);
   DockTK::nr::jacobi(c_nrc.rowPointers(),4,d.dataFirst()-1,v_nrc.rowPointers(),&nrot); //jacoby diagonalization
   DockTK::nr::eigsrt(d.dataFirst()-1,v_nrc.rowPointers(),4); //sort eigenvalues

   //extract the desired quaternion
   vect4 q;
   /* Ponder's jacobi sorts eigenvalues in the order opposite to the "eigsrt" routine from NRC.
   We switch the order.
   Ponder:
   q(1) = v(1,4);
   q(2) = v(2,4);
   q(3) = v(3,4);
   q(4) = v(4,4);
   Eigsrt:
   */
   q(0) = v(0,0);
   q(1) = v(1,0);
   q(2) = v(2,0);
   q(3) = v(3,0);

   //DEBUG:
   //cout << "d=" << d << endl;
   
   //assemble the rotation matrix that superimposes molecules
   //matrix is transposed compared to Fortran version
   typename QuatFitTraits<T_num>::Matr3 rot;
   //using namespace ::DockTK::math;
   rot(0,0) = DockTK::Math::pow2(q(0)) + DockTK::Math::pow2(q(1)) - DockTK::Math::pow2(q(2)) - DockTK::Math::pow2(q(3));
   rot(0,1) = 2.00 * (q(1) * q(2) - q(0) * q(3));
   rot(0,2) = 2.00 * (q(1) * q(3) + q(0) * q(2));
   rot(1,0) = 2.00 * (q(2) * q(1) + q(0) * q(3));
   rot(1,1) = DockTK::Math::pow2(q(0)) - DockTK::Math::pow2(q(1)) + DockTK::Math::pow2(q(2)) - DockTK::Math::pow2(q(3));
   rot(1,2) = 2.00 * (q(2) * q(3) - q(0) * q(1));
   rot(2,0) = 2.00 * (q(3) * q(1) - q(0) * q(2));
   rot(2,1) = 2.00 * (q(3) * q(2) + q(0) * q(1));
   rot(2,2) = DockTK::Math::pow2(q(0)) - DockTK::Math::pow2(q(1)) - DockTK::Math::pow2(q(2)) + DockTK::Math::pow2(q(3));

   return rot;
}

} } } // namespace DockTK { namespace Geom::Align
