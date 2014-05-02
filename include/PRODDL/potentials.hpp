//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PRODDL_POTENTIALS_H__
#define AT_PRODDL_POTENTIALS_H__

#include "PRODDL/Common/math.hpp"

#include "PRODDL/types.hpp"

#include "PRODDL/Geom/pairdist.hpp"

#include "PRODDL/Common/g_options.hpp"

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/logger.hpp"

#include <vector>

#ifdef PRODDL_CONTRIB
#   include "PRODDL/Contrib/potential_step.hpp"
#endif


// "Parametrized typedef" class that wraps classes for intercation potentials

namespace PRODDL { 

  template<typename T_num>
  class
  Potentials {

  public:


    typedef typename Types<T_num>::Point Point;
    typedef typename Types<T_num>::Points Points;
    typedef typename Types<T_num>::Matrix Matrix;
    typedef typename Types<T_num>::IntMatrix IntMatrix;
    typedef typename Types<T_num>::Ints Ints;
    typedef typename Types<T_num>::Uints Uints;
    typedef T_num Float;
    typedef typename Types<T_num>::Floats Floats;
    typedef Geom::Points::Neighbors<T_num,3> NeighborsType;
    typedef typename NeighborsType::PartitionedPoints2 PartPoints;

    // Faster implementaion of PartPoints
    typedef typename Geom::Points::PartitionedPoints<T_num,int,3> PartPointsFast;
    typedef Geom::Points::PositionExtractorByArrayIndex<T_num,3> PosExtractor;
    typedef Geom::Points::SphereIter<PartPointsFast,PosExtractor> SpIter;

    typedef typename NeighborsType::VIPair VIPair;
    typedef typename NeighborsType::IPair IPair;
    typedef typename NeighborsType::fvect fvect;



#   include "PRODDL/potentials_inc.hpp"

#   include "PRODDL/potentials_types.hpp"

#   include "PRODDL/potentials_total.hpp"

#   include "PRODDL/potentials_total_rot.hpp"


  }; // class Potentials { 

} // namespace PRODDL

#endif // AT_PRODDL_POTENTIALS_H__
