//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef GRID_PROJECT_H__
#define GRID_PROJECT_H__


#include "PRODDL/Grid/grid.hpp"

#include "PRODDL/Common/nd_index_iter.hpp"

#include "PRODDL/Common/bz_vect_ext.hpp"

#include "PRODDL/Common/common_types.hpp"

namespace PRODDL { namespace Grid {


  // Class Projector is used as a 'parametrized namespace'.
  // The purpose of this class is to provide common typedefs and simplify declarations of member
  // function pointers (for Python interface, for instance)
  // Template parameters:
  // T_val - value type of the grid
  // T_coord - coordinate type of the grid
  // n_dim - dimentionality of the grid and space

  template<typename T_val,typename T_coord,int n_dim> 
  class Projector {

  public:

    typedef Grid<n_dim,T_val,T_coord> Grid_t;
    typedef typename Grid_t::SpatialCoord SpatialCoord;
    typedef typename Grid_t::LogicalCoord LogicalCoord;

    typedef typename common_types::num_vector_type<SpatialCoord>::Type SpatialCoords;

    // The helper class that provides operations common to several types of projection
    // functions (different projection functions project value of radial potential,
    // gradient of radial potential etc)
    // This class defines the rectangular subdomain region corresponding to 'cutoffRadius' sphere
    // around 'center' of potential, and allows to iterate through all cells of that region.
    // The object of this class is expected to be created on the stack within
    // the function that does the projection.

    class GridRegionProjectionHelper {
    public:

      typedef PRODDL::index_mover<n_dim> ind_mover;
      enum { lBound = 0, uBound = 1 };

      // cutoffRadius^2
      T_coord cutoffRadiusP2;

      // region around 'center'
      typename Grid_t::LogicalDomain region;

      // current cell index during iteration
      typename Grid_t::LogicalCoord indCurr;
    
      // ctor
      // Arguments:
      // 'center' - potential center
      // 'cutoffRadius' - cutoff radius for the projection.
      // 'grid' - this grid object will keep the projection.
      // Results:
      // After this ctor, next() member can be called in a loop untill it returns false. That will advance
      // 'indCurr' cell index. Initial value of 'indCurr' is before the first valid element.

      GridRegionProjectionHelper(const typename Grid_t::SpatialCoord& center,
				 T_coord cutoffRadius,
				 Grid_t& grid) :

	cutoffRadiusP2(cutoffRadius * cutoffRadius) 
      {
      
	// TODO: (if necessary for performance reasons)
	// The following method of calculating rectangular region will include some cells which 
	// have centers outside of cutoffRadius (e.g. when 'center' is to the left of its cell's center,
	// inclusion of the rightmost cells will be excessive). This can be changed to the more
	// elaborate scheme if the projection stage becomes time critical part of overall calculation
	// (which is unlikely).

	region(lBound) = grid.getGeometry().toLogical(center - cutoffRadius);
	region(uBound) = grid.getGeometry().toLogical(center + cutoffRadius);

	// potentially crop the region to make sure all indices within it are for valid grid cells
      
	grid.cropSubDomain(region);

	region(uBound) += 1; // index_mover needs open range

	ind_mover::before_first(indCurr,region(lBound));
      }

      // move to the next cell in 'region'
    
      inline bool next() {
	return ind_mover::next(indCurr,region(lBound),region(uBound));  
      }
      
    }; // class GridRegionProjectionHelper


    // Additively project radial 'field' on a 'grid' within the 'field.getRadius()' around field's 'center' (source). Only cells
    // whose centers are within 'field.getRadius()' from fields's source 'center' are affected. 
    // RadialField object must have: 
    // - a method f2(T_coord rP2) that takes SQUARE of a distance as argument and returns the value of the field at this 
    // distance from the source of the field. In each call to 'field', distance argument is guaranteed to be within 'field.getRadius()' 
    // center, so RadialField can be optimized relying on that condition (if the field is a step function, for instance).
    // - a method getRadius() that returns the distance from the field's center within which the field has non-zero values. The area
    // outside this radius is considered to be zero and ignored during projection.
    // Template parameters:
    // RadialField - potential type (see above for requirements)
    // Arguments:
    // 'center' - field's source (e.g. point charge location for electrostatic potential)
    // Results:
    // Values in grid cells within 'field.getRadius()' from 'center' are INCREMENTED by the values of the field at their centers.
    // Complexity: (field.getRadius()/grid_step)^n_dim
  
    template<class RadialField> 
    class ProjectorRadialField {

    public:

      typedef typename common_types::num_vector_type<RadialField>::Type RadialFields;

      ProjectorRadialField() :
	m_pgrid(0),
	m_minDist2(0)
      {}

      void init(Grid_t& grid,T_coord minDist) {
	m_pgrid = &grid;
	m_minDist2 = minDist * minDist;
      }

      void projectField(RadialField& field, 
			const SpatialCoord& center)
      {
	typedef GridRegionProjectionHelper _ProjectionHelper;
	Grid_t& grid = *m_pgrid;
	_ProjectionHelper projectionHelper(center,field.getRadius(),grid);
	while( projectionHelper.next() ) {
	  SpatialCoord posCurr = grid.getGeometry().toSpatial(projectionHelper.indCurr);
	  T_coord rP2 = blitz_ext::dotSelf(posCurr-center);
	  if( rP2 <= projectionHelper.cutoffRadiusP2 ) {
	    if( rP2 < m_minDist2 ) {
	      rP2 = m_minDist2;
	    }
	    grid.getGridArray()(projectionHelper.indCurr) += field.f2(rP2);
	  }
	}
      }

      void projectFields(RadialFields& fields, 
			 const SpatialCoords& centers,
			 bool zeroGrid = true) {

	Grid_t& grid = *m_pgrid;
	if( zeroGrid ) {
	  grid = 0;
	}

	for( int i = 0; i < fields.size(); i++ ) {

	  projectField(fields(i),centers(i));

	}

      }

    protected:

      Grid_t * m_pgrid;
      T_coord m_minDist2;

    }; // class ProjectorRadialField

  }; // class Projector

} } // namespace PRODDL::Grid


#endif // GRID_PROJECT_H__
