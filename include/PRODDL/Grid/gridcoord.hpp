//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef __ATGRIDCOORD_H
#define __ATGRIDCOORD_H


//Conversion to/from spacial/index coordinates
//Made with small modifications from Blitz geometry.h.
//Currently only toLogical() method added. We don't
//inherit from original Blitz class because Blitz/geometry.h
//is undocumented.

#include <limits> //for blitz/funcs which use std:: but include limits.h
#include <blitz/tinyvec2.h>
//#include <blitz/applics.h> //for pow()
#include <blitz/funcs.h> //for pow()

namespace PRODDL { namespace Grid {

  namespace Geometry {
    using blitz::TinyVector;
    using blitz::pow2;
    using blitz::pow3;
    
    typedef double T_defaultSpatialCoordinate;
    
    template<int N_dim, class T = T_defaultSpatialCoordinate>
      class UniformOrthoGeometry {
      public:
    };
    
    //Note: If we use this object to allocate the grid array:
    //UniformCubicGeometry geom(step,space_lower_bound,ind_offset)
    //grid = Array(ind_offset,geom.toLogical(space_upper_bound)),
    //and use same object later to obtain indices into array:
    //ind = geom.toLogical(space_point); grid(ind), then
    //space_point is allowed to be in 
    //(space_lower_bound, space_upper_bound]. Left limit is exact,
    //while right limit can be up to space_upper_bound+h/2 if step is not
    //rounf divisor of space interval.
    // The method used to convert space coords into logical coords allows
    // for arrays with negative offsets. That achieved by small performance
    // trade-off nessesary to use int() conversion instead of slower
    // floor() function. See toLogical() and setup() for details.
    //Bins are: 
    //(space_lower_bound, space_lower_bound + step)
    //[space_lower_bound + step, space_lower_bound + step*2)
    //...
    //[space_lower_bound + step*(n - 1), space_lower_bound + step*n].
    //toSpatial() returns the center of each bin.
    //To enumerate bins' centers, use something like:
    //for(x = geom.toSpatial(grid.lbound()); x <= geom.toSpatial(grid.ubound()); x += step)
    //To get domain of spatial coords, use SpatialDomainLower() and SpatialDomainUpper().
    //They will return the left boundary of the first bin and the right boundary of the
    //last bin correspondingly.

template<int N_dim, class T = T_defaultSpatialCoordinate>
  class UniformCubicGeometry {
  //spatial step
  TinyVector<T,N_dim> h_;
  //precomputed derivatives of h_
  TinyVector<T,N_dim> recip_h_;
  TinyVector<T,N_dim> recip2_h_;
  TinyVector<T,N_dim> recip3_h_;
  //the center of first cell
  TinyVector<T,N_dim> zero_;
  //coord base used to convert into indices
  TinyVector<T,N_dim> zero_start_;
  //coord base used to convert into indices in toLogical()
  TinyVector<T,N_dim> zero_h_;
  //coord base used to convert indices into cell centers
  TinyVector<T,N_dim> zero_start_center_;
  //indices will be made having this offset
  TinyVector<int,N_dim> offset_;
  //used in toLogical()
  TinyVector<int,N_dim> offset_1_;
    
  public:
    typedef T T_coord;

    typedef UniformCubicGeometry<N_dim,T_coord> Self;
    
    UniformCubicGeometry()
      {
        h_ = 0.0;
        recip_h_ = 0.0;
        recip2_h_ = 0.0;
        recip3_h_ = 0.0;
	offset_ = 0;
	offset_1_ = 0;
	zero_start_center_ = 0;
	zero_start_ = 0;
	zero_ = 0;
	zero_h_ = 0;
      }
    
    
    UniformCubicGeometry(  TinyVector<T,N_dim> spatialStep, 
			   TinyVector<T,N_dim> zeroCoordinates,
			   TinyVector<int,N_dim> indOffset  )
      {   
        h_ = spatialStep;
        zero_ = zeroCoordinates;
	offset_ = indOffset;
        setup();
      }    

    UniformCubicGeometry(  TinyVector<T,N_dim> spatialStep, 
			   TinyVector<T,N_dim> zeroCoordinates)
      {   
        h_ = spatialStep;
        zero_ = zeroCoordinates;
	offset_ = 0;
        setup();
      }    


    // Geometry will be defined as a rectangular region with a center 
    // approximately at 0, and a logical size 'logicalSize'

    UniformCubicGeometry(  TinyVector<T,N_dim> spatialStep, 
			   TinyVector<int,N_dim> logicalSize)
      {   
        h_ = spatialStep;
	TinyVector<int,N_dim> logical_zero;
	logical_zero = 0;
	TinyVector<T,N_dim> spatial_zero;
	spatial_zero = 0;
	Self geom(h_,spatial_zero);
	TinyVector<T,N_dim> s_lBound = geom.toSpatialLeft(logical_zero);
	TinyVector<T,N_dim> s_uBound = geom.toSpatialRight(logicalSize-1);
	TinyVector<T,N_dim> sSize = s_uBound - s_lBound;
        zero_ = -1. * sSize/2;
	offset_ = 0;
        setup();
      }    

    
    TinyVector<T,N_dim> toSpatial(TinyVector<int,N_dim> logicalCoord) const
      {
        return zero_start_center_ + h_ * logicalCoord;
      }

    TinyVector<T,N_dim> toSpatialDiff(TinyVector<int,N_dim> logicalDiff) const
      {
        return h_ * logicalDiff;
      }
    
    const TinyVector<T,N_dim>& spatialStep() const
      { return h_; }
    
    const TinyVector<T,N_dim>& recipSpatialStep() const
      { return recip_h_; }

    const TinyVector<T,N_dim>& recipSpatialStepPow2() const
      { return recip2_h_; }
    
    
private:
    void setup()
      {
        recip_h_ = 1.0 / h_;
        recip2_h_ = 1.0 / h_ / h_;
        recip3_h_ = 1.0 / h_ / h_ / h_;
	zero_start_ = zero_ - h_*offset_;
	zero_h_ = zero_ - h_;
	zero_start_center_ = zero_start_ + h_/2.0;
	offset_1_ = offset_ - 1;
      }
    
  public:
    
    //rounds given spatial coord to the corresponding grid point. For anything less
    //than spatialZero() it will return something less than logicalZero(), but the exact
    //value does not hold much sence in that case and can be used only
    //as a flag. Use slower toLogicalUnlimited() after you check the return value
    //from toLogical() to get the correct index for coords below spatialZero().

    TinyVector<int,N_dim> toLogical(TinyVector<T,N_dim> spatialCoord) const
    {
      // the expession below garantees that arguments between spatialZero()-1 and spatialZero()
      // will be rounded to logicalZero()-1 rather than to logicalZero(). This trick
      // will make such coords to be 'out-of-bound', as they should be.

      return TinyVector<int,N_dim>((spatialCoord - zero_h_)/h_) + offset_1_;
    }

    TinyVector<int,N_dim> toLogicalDiff(TinyVector<T,N_dim> spatialDiff) const
    {
      return TinyVector<int,N_dim>(spatialDiff/h_);
    }

    //returns correct results for any coords, including less than zero()
    TinyVector<int,N_dim> toLogicalUnlimited(TinyVector<T,N_dim> spatialCoord) const
    {
      TinyVector<T,N_dim> res((spatialCoord - zero_start_)/h_);
      for(int i=0; i < N_dim; i++) {
          res[i] = blitz::floor(res[i]);
      }
      return res;
    }
    
    const TinyVector<T,N_dim>& spatialZero() const
      { return zero_; }

    const TinyVector<int,N_dim>& logicalZero() const
      { return offset_; }
    
    //return the leftmost spatial coordinate of the bin corresponding to the
    //given logical coordinate. This function has unlimited domain for it's argument
    //(i.e., it can be negative)
    TinyVector<T,N_dim> toSpatialLeft(TinyVector<int,N_dim> logicalCoord) const
      {
	return zero_start_ + h_ * logicalCoord;
      }

    //return the rightmost spatial coordinate (exclusive) of the bin corresponding to the
    //given logical coordinate. This function has unlimited domain for it's argument
    //(i.e., it can be negative)
    TinyVector<T,N_dim> toSpatialRight(TinyVector<int,N_dim> logicalCoord) const
      {
	return zero_start_ + h_ + h_ * logicalCoord;
      }

};
 
 template<int N_dim, class T = T_defaultSpatialCoordinate>
   class TensorProductGeometry {
   public:
 };
  } //namespace Geometry
} } //namespace PRODDL::Grid
  
#endif // __ATGRIDCOORD_H
  
