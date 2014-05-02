//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_GEOM_QUATFIT_H__
#define AT_GEOM_QUATFIT_H__

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

#include <blitz/tinymat.h>

namespace PRODDL { namespace Geom { namespace Align {

template<typename T_num> struct QuatFitTraits {
  typedef blitz::TinyMatrix<T_num,3,3> Matr3;
};

template<typename T_num,class PointsInputIterator1,class PointsInputIterator2,class WeightsInputIterator>
typename QuatFitTraits<T_num>::Matr3 
quatFit(PointsInputIterator1 points1_first,
	PointsInputIterator1 points1_last,
	PointsInputIterator2 points2_first,
	WeightsInputIterator weights_first);


} } } // namespace PRODDL { namespace Geom::Align

#include "PRODDL/Geom/quatfit.cc"

#endif // AT_GEOM_QUATFIT_H__

