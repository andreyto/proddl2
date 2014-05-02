//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_GEOM_POINTS_H__
#define PRODDL_GEOM_POINTS_H__

#include "PRODDL/Common/math.hpp"

#include "PRODDL/types.hpp"

#include "PRODDL/Geom/pairdist.hpp"

#include "PRODDL/Geom/bounding.hpp"

#include "PRODDL/Common/g_options.hpp"

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/logger.hpp"

#include <vector>



// "Parametrized typedef" class that wraps classes for space points manipulation

namespace PRODDL { 

  template<typename T_num>
  class
  Points {

  public:


    typedef typename Types<T_num>::Point Point;
    typedef typename Types<T_num>::Points Points;
    typedef typename Types<T_num>::Ints Ints;
    typedef T_num Float;
    typedef typename Types<T_num>::Floats Floats;
    typedef Geom::Points::Neighbors<T_num,3> NeighborsType;
    typedef typename NeighborsType::PartitionedPoints2 PartPoints2;
    typedef PartitionedPoints<T_num,int,3> PartPointsInt;
    typedef PositionExtractorByArrayIndex<T_num,3> PosExtractorInt;
    typedef SphereIter<PartPointsInt,PosExtractorInt> SphereIterInt;
    typedef typename PRODDL::common_types::num_vector_type<Points> PointsSet;
    typedef typename NeighborsType::VIPair VIPair;
    typedef typename NeighborsType::IPair IPair;
    typedef typename NeighborsType::fvect fvect;

    typedef Geom::Bounding::Box<T_num> BoundingBox;
    typedef Geom::Bounding::Diameter<T_num> BoundingDiameter;
    typedef typename BoundingBox::PointPair PointPair;


    // NOT USED CURRENTLY FOR ANYTHING, NOTHING TO INCLUDE HERE

  }; // class Points { 

} // namespace PRODDL

#endif // PRODDL_GEOM_POINTS_H__
