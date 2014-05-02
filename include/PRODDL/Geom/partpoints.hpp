//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef  AT_PARTPOINTS_H__
#define  AT_PARTPOINTS_H__

/* Class to organize efficient loop over pairs of points
   located within predefined cutoff distance from each other by
   assigning them to cells on a grid. For a given point, iterator-like object
   can be obtained and used in a usual for(...) syntax to loop over all points
   in the nearest cells. Dimensionality of a space is a template parameter.

   Update: See also examples with the more compact iterator expressions
   using SphereIter class further in this file.

   //////////
   // Example 1:
   // 'aPoints' is array of Point 3D points distributed within region [(-100,-100,-100),(100,100,100)].
   // The task is to find pairs of points from 'aPoints' located within cutoff distance 1.7 from each
   // other.
   const double cutoff = 1.7;
   const double cutoffP2 = cutoff*cutoff;
   typedef PartitionedPoints<double,int,3> PartPoints;
   PartPoints partPoints(Point(-100,-100,-100),Point(100,100,100),cutoff);
   for(int i_point = 0; i < aPoints.size(); i_point++) {
   Point& point = aPoints[i_point];
   // get CellIndex object for Point 'point' - CellIndex references
   // the cell to which 'point' would belong
   PartPoints::CellIndex cellInd = partPoints.getCellIndex(point);
   // iterate over all points in the neghbor cells. All points that might be within cutoff distance
   // from 'point' are garanteed to be included in this set
   for(PartPoints::SubDomainIter iterNeighb = cellInd.getNeighbors(); iterNeighb.not_end(); iterNeighb.next()) {
   int i_point_other = *iterNeighb;
   Point& point_other = aPoints[i_point_other];
   if( bz_ext::dotSelf(point - point_other) < cutoffP2 ) {
   cout << "Points (" << i_point << ", " << i_point_other << ") are within cutoff distance" << endl;
   }
   }
   // now store this point's index in corresponding grid cell
   cellInd.insert(i_point);
   }
   // Note that each point is inserted into corresponding grid cell after looping over already stored points from
   // neighboring cells. This ensures that each pair is counted only once.
   // End Example 1.

   ///////////

   // Example 2:
   // 'aPoints1' and 'aPoints2' are two arrays of Point 3D vectors distributed within region 
   // [(-100,-100,-100),(100,100,100)].
   // The task is to find pairs of points (i,j) such that 'i' belongs to  'aPoints1', 'j' belongs to 'aPoints2' 
   // and distance (i,j) is less than cutoff distance 1.7.
   const double cutoff = 1.7;
   const double cutoffP2 = cutoff*cutoff;
   typedef PartitionedPoints<double,int,3> PartPoints;
   PartPoints partPoints(Point(-100,-100,-100),Point(100,100,100),cutoff);
   // first, insert all points from set 2 into the grid
   for(int i_point = 0; i < aPoints2.size(); i_point++) {
   Point& point = aPoints2[i_point];
   // store into 'cellInd' CellIndex for Point 'point' - CellIndex references
   // the cell to which 'point' would belong
   PartPoints::CellIndex cellInd = partPoints.getCellIndex(point);
   cellInd.insert(i_point);
   }
   // Now, for each point from set 1 find neighbor points from the grid.
   for(int i_point = 0; i < aPoints1.size(); i_point++) {
   Point& point = aPoints1[i_point];
   // store into 'cellInd' CellIndex for Point 'point' - CellIndex references
   // the cell to which 'point' would belong
   PartPoints::CellIndex cellInd = partPoints.getCellIndex(point);
   // iterate over all points in the neghbor cells. All points that might be within cutoff distance
   // from 'point' are garanteed to be included in this set
   for(PartPoints::SubDomainIter iterNeighb = cellInd.getNeighbors(); iterNeighb.not_end(); iterNeighb.next()) {
   int i_point_other = *iterNeighb;
   Point& point_other = aPoints2[i_point_other];
   if( bz_ext::dotSelf(point - point_other) < cutoffP2 ) {
   cout << "Points (" << i_point << ", " << i_point_other << ") are within cutoff distance" << endl;
   }
   }
   }
   // End Example 2.

   ///////////
   */

//#include <list>
// vector is slightly more efficient for small
// datatypes (typical case)
#include <vector>
#include <algorithm>
#include "PRODDL/Common/bz_vect_ext.hpp"

#include "PRODDL/Common/nd_index_iter.hpp"

#include "PRODDL/Grid/grid.hpp"

#include "PRODDL/Blitz/bzalgor.hpp"

namespace PRODDL { namespace Geom { namespace Points {

  // Template parameters:
  // _T_num - numeric type to use (e.g. double).
  // _PointData - objects of this type will be stored in a lists assigned to grid cells
  // (e.g. Point - 3d vector of position, or int - as index to some array of 3d vectors).
  // _n_dim - points are defined in _n_dim-space.

  template<typename _T_num,class _PointData,int _n_dim>
  class PartitionedPoints {

    friend class SubDomainIter;

  public:

    enum { N_dim = _n_dim };
    typedef _T_num T_num;
    typedef _PointData PointData;

    // type of lists of points stored at grid cells
    typedef std::vector<PointData> Cell;

    // type of grid
    typedef ::PRODDL::Grid::Grid<_n_dim,Cell,T_num> GridCut;

    // N_dim points in space
    typedef typename GridCut::SpatialCoord Point;

    // Index of grid elements
    typedef typename GridCut::LogicalCoord Index;

    // Nested class SubDomainIter provides the methods to organize the loop
    // over points located inside a NxN-cells grid area
    // centered at 'indCenter' cell where N = 2*radius+1 (default radius is 1),
    // or over grid area defined by logical domain (see descriptions of ctors).
    // The objective is to use a single for(...) loop to iterate over all points
    // in neighbor cells.
    // The class methods must be used in specific order. Example:
    // for(SubDomainIter iter(partitionedPoints,indCenter); iter.not_end(); iter.next()) {
    //   if(dot(*iter,pointCenter) < cutOffP2) {
    //      .....
    //   }
    // }
    // End Example.
    // 'iter' is dereferenceable only after call to iter.not_end(). iter.next() can be called only
    // once after each call to iter.not_end().
    // Method init(indCenter) can be used instead of constructor in the example above if already
    // created (and, possibly, used) object is available.
    // TODO: add missing components to support full std::iterator concept: iterator_category, const-nonconst,
    // operator==, postfix and prefix++ etc. Comparison will cost one more 'if' statement (curr != end);
    // postfix++ will cost returning a copy of self. Const-nonconst can be handled by atalgor.h:same_const<>
    // mechanism. Implementation of Std::deque<> can be used as an example. Note: this is a lot of work, deffer it untill
    // we really need it somewhere.

    class SubDomainIter {

    public:
      typedef typename Cell::iterator cell_iterator;
      typedef typename std::iterator_traits<cell_iterator>::reference reference;
      typedef typename std::iterator_traits<cell_iterator>::pointer pointer;

    protected:
      typedef ::PRODDL::index_mover<N_dim> ind_mover;

    public:
      

      // default ctor for constructing arrays of these objects
      SubDomainIter():
	pOwner(0)
      { }

      // create iterator object for cells' cube with side=2*radius+1 centered at indCenter cell

      SubDomainIter(PartitionedPoints& owner, const Index& indCenter,const int radius=1): 
	pOwner(&owner) {
	  init(indCenter,radius);
	}

      // create iterator object for cells with indices in [_indLow,_indUp] (note the closed range)

      SubDomainIter(PartitionedPoints& owner, const Index& _indLow,const Index& _indUp): 
	pOwner(&owner), indLow(_indLow), indUp(_indUp) {
	  ++indUp; // make open range expected by 'ind_mover'
	  reset();
	}

      // create object for the entire grid if 'doInit'
      // or just assign 'owner' and leave the rest
      // as it is until some 'init()' method is called
      // explicitly

      SubDomainIter(PartitionedPoints& owner, bool doInit = true) :
	pOwner(&owner) {
	  if( doInit ) {
	    indLow = pOwner->_getGrid().getLogicalDomain()(0);
	    indUp = pOwner->_getGrid().getLogicalDomain()(1);
	    indUp += 1; // make open range expected by 'ind_mover'	
	    reset();
	  }
	}

      // see corresponding ctor
      void init(const Index& indCenter,const int radius=1) {
	indLow = indCenter - radius;
	indUp  = indCenter + radius;
	bool is_cropped = false;
	// the following sequence will
	// crop index limits indLow and indUp to the corresponding nearest cells if 
	// the initial values are out of bounds:
	GridCut& grid = pOwner->_getGrid();
	blitz_ext::crop(indLow,grid.getLogicalDomain()(0),grid.getLogicalDomain()(1),is_cropped);
	is_cropped = false;
	blitz_ext::crop(indUp,grid.getLogicalDomain()(0),grid.getLogicalDomain()(1),is_cropped);
	indUp+=1; // make open range expected by 'ind_mover'
	// end of cropping indLow and indUp
	reset();
      }

      // see corresponding ctor
      void init(const Index& _indLow, const Index& _indUp) {
	indLow = _indLow;
	indUp = _indUp;
	indUp+=1; // make open range expected by 'ind_mover'
	reset();
      }

      // reset iterator to the initial position defined by indLow and indUp values, which
      // must be already set up before calling this function. Must be: [indLow,indUp)

      void reset() {
	ind_mover::before_first(indCurr,indLow);
	// set initial values of cell iterators such that iterCell == iterCellEnd and both are valid
	// iterators in a sense the iterCell != iterCellEnd is a valid operator (where they point at
	// this moment does not matter). That is essential for the correct first call to 'not_end()'
	// method implementation
	Cell& cell = pOwner->grid.getGridArray()(indLow);
	iterCellEnd = cell.end();
	iterCell = iterCellEnd;
      }

      // check if current cell iterator is within the valid range for the current cell,
      // otherwise advance to the next non-empty cell. If this function returned true,
      // current iterator always points to a valid PointData object. After this function
      // returned false, it should not be called anymore untill 'init()' is called.

      bool not_end() {
	// this condition is true most of the time
	// for the non-empty cell
	if( iterCell != iterCellEnd ) {
	  return true;
	}
	// find the next non-empty cell and initialize cell iterators
	while( ind_mover::next(indCurr,indLow,indUp) ) {
	  Cell& cell = pOwner->grid.getGridArray()(indCurr);
	  iterCellEnd = cell.end();
	  iterCell    = cell.begin();
	  if(iterCell != iterCellEnd)
	    return true;
	}
	return false;
      }

      // This function can be called once only after call to 'not_end()' returned true.

      void next() {
	++iterCell;
      }

      // dereferencing this iterator will return a reference to current PointData object

      reference operator*() const {
	return *iterCell;
      }

      pointer operator->() const { return &(*iterCell); }

      PartitionedPoints& getOwner() const {
	return *pOwner;
      }

    protected:

      Index indLow, indUp, indCurr;
      PartitionedPoints *pOwner;
      cell_iterator iterCell,iterCellEnd;

    }; // end of class SubDomainIter


  public:

    // small proxy object to provide operations typically
    // performed while indexing PartitionedPoints object.
    // This object must be produced by PartitionedPoints indexing operation. 

    class CellIndex {

    public:

      CellIndex(): _p_cell(0), _pOwner(0) {
      }

      CellIndex(PartitionedPoints& owner): _p_cell(0), _pOwner(&owner) {
      }

      Cell& getCell() { return *_p_cell; }

      const Index& getIndex() const { return _ind_cell; }

      // insert PointData into this cell

      void insert(const PointData& x) {
	_p_cell->push_back(x);
      }

      // obtain iterator to loop over PointData objects stored in
      // this cell and its neighbors. All points that might be within cutoff distance
      // from any point assigned to this cell are guaranteed to be included in this set

      SubDomainIter getNeighbors() {
	return SubDomainIter(*_pOwner,_ind_cell);
      }

    public:
      
      // this members are made public to let PartitionedPoints access
      // them directly when creating instances of this class. They
      // are implementation specific and should not be accessed by user code.

      Cell *_p_cell;
      Index _ind_cell;
      PartitionedPoints *_pOwner;
    };


  public:

    // Default ctor creates empty object, which must be intialized later by
    // a call to init(...) method.

    PartitionedPoints() {}

    // Construct object based on lower and upper corners of a rectangular space region, which
    // the points are expected to occupy [lBoundS,uBoundS] and a cutoff value defining the
    // neigbors 'cutoff'. Internally, the grid will be created with spatial step of 'cutoff' and 
    // PointData objects supplied with their Point coordinates will be stored in corresponding
    // grid cells. Therefore, for a given point, finding all its neighbors takes time O(n_cell) where
    // n_cell is a number of points in a cell. The choice of a spatial domain must correspond to the properties
    // of the actual points distribution to provide the best performance and memory usage.
    // Ideal case is when the points are distributed uniformly within the supplied rectangular domain
    // [lBoundS, uBoundS]. Then time for N points is O(N * N / N_cells) where N_cells is a number of cells in
    // a grid. For instance, if the task is to find VDW contacts in a globular protein, cutoff is ~1.7 and N_cells is
    // O(N), which gives time O(N).
    // Worst case scenario is when all the points are located within a single grid cell. When performance
    // is quadratic over a number of points like in a simple double loop algorithm plus some amortized constant
    // time penalties for the greater code complexity.
    // It is allowed for a given point to fall outside of the supplied grid domain. 
    // Such point will be assigned to the nearest border cell. Thus, even
    // if the points are distributed in much larger area than the supplied domain but in such manner that
    // they are assigned uniformly to different border cells (e.g. points are uniformly distributed over a larger area),
    // then there will be still a significant performance gain (O(N_border_cells)).

    PartitionedPoints(const Point& lBoundS,const Point& uBoundS, T_num cutoff) {
      init(lBoundS, uBoundS, cutoff);
    }

    // reinitialize the object - all content is erased, grid memory is reallocated. 
    // For arguments see ctor description.
    void init(Point lBoundS, Point uBoundS, T_num cutoff) {
      // we check if uBoundS-lBoundS will lead to co-planar grid vertices and extend
      // their range by 'cutoff' to remove such condition
      for(int i = 0; i < N_dim; i++) {
	if(uBoundS(i) - lBoundS(i) <= cutoff) {
	  uBoundS(i) += cutoff/2.;
	  lBoundS(i) -= cutoff/2.;
	}
      }
      grid.init(lBoundS,uBoundS,cutoff,false);
      grid.setDomainAction(GridCut::domain_crop);
    }


    // Clear all data about previously accumulated points.
    // Grid structure is not changed.
    // After this method is called, the new loop of inserting points can be started
    // immediately.
    void zapPoints() {
      grid.getGridArray() = Cell();
    }

    T_num getCutoff() const {
      return grid.getSpatialStep()[0];
    }

    T_num getCutoffP2() const {
      T_num cutoff = getCutoff();
      return cutoff*cutoff;
    }


    // Return the CellIndex corresponding to the cell into which
    // Point 'point' is projected. The returned object can
    // be used to access cell content and to iterate over
    // all candidates for 'point''s neighbors (through 'getNeighbors()' method).

    CellIndex getCellIndex(const Point& point) {
      CellIndex cell_ind;
      getCellIndex(point,cell_ind);
      return cell_ind;
    }

    // Same as above but reinitializes existing CellIndex object.

    void getCellIndex(const Point& point,CellIndex& cell_ind) {
      cell_ind._pOwner = this;
      grid.indAt(point,cell_ind._ind_cell);
      cell_ind._p_cell = &(grid.getGridArray()(cell_ind._ind_cell));
    }

    //insert PointData object with a given position

    void insert(const Point& pos, const PointData& x) {
      grid.at(pos).push_back(x);
    }
    

  public:

    // The methods in this section made public for the sake
    // of closely related classes like SphereIter.
    // They should be used by the user code.

    // return internal grid object
    GridCut& _getGrid() {
      return grid;
    }

  protected:

    // PointData objects are stored here
    GridCut grid;

  }; // class PartitionedPoints




  // SphereIter class is designed to provide syntactically
  // simple code when one needs to iterate over all
  // PointData objects within a given sphere around a given
  // point in space.
  // Because SphereIter needs to compute the distance between the center
  // point and PointData objects stored inside the PartitionedPoints object,
  // we have to provide a generic way to compute this distance
  // (position might not be stored inside the PointData object -
  // PointData might be just an index into some array).
  // Thus, SphereIter is a template class whose ctor accepts
  // a function object TPositionExtractor pos, that when called
  // as pos(const PointData& pd) must return const Point& with
  // the position of 'pd'.
  // Because SphereIter is not a small object, it must be first
  // created, and then reinitialized for each new 'pointCenter'
  // with a call to init() method.
  // As with the SubDomainIter class, the sequence of calls to
  // the loop-related methods must be strictly maintained, as in
  // the example below (next() must be called when and only when
  // not_end() returned truth; calling not_end() w/o the interleaving
  // call to next() is an error).
  // Usage Example:
  // struct MyPositionExtractor {
  //     const Point& operator()(const PointData& o) const {
  //         return o.pos;
  //     }
  // };
  // ... ...
  // //If radius is <0 or omitted, cutoff of partitionedPoints
  // //will be used.
  // SphereIter<MyPositionExtractor> iter(partitionedPoints,
  //                                      MyPositionExtractor(),
  //                                      radius);
  // while(pointCenter = ...) {
  // for(iter.init(pointCenter); iter.not_end(); iter.next()) {
  //
  //     // o is always within 'radius' distance from 'pointCenter'
  //     const PointData& o = *iter;
  //     T_num r2 = iter.r2();
  //     ... Do something with o and r2 ...
  //
  // }
  // //optionaly:
  // iter.insert(PointData(pointCenter,whatever_else));
  // } // end while
  // End Example.
  // Alternative call sequence Example:
  // while(pointCenter = ...) {
  // for(iter.init(pointCenter); iter.not_end_lazy(); iter.next()) {
  //
  //     // o is not nesessarily within 'radius' distance from 'pointCenter',
  //     // the distance calculation is delayed.
  //     const PointData& o = *iter;
  //     // check for some other (faster than distance calculation) conditions
  //     // too see if this 'o' element should be considered (e.g. pairwise
  //     // interaction exclusion matrix), and only then check/compute the distance:
  //     if ( useThisElement(o) && iter.check_dist() ) {
  //     {
  //         T_num r2 = iter.r2();
  //         ... Do something with o and r2 ...
  //     }
  //
  // }
  // } // end while
  // End Example.



  template<class TPartitionedPoints, class TPositionExtractor>
  class SphereIter {

  public:
      
    typedef typename TPartitionedPoints::T_num T_num;
    typedef typename TPartitionedPoints::PointData PointData;
    typedef typename TPartitionedPoints::Point Point;
    typedef typename TPartitionedPoints::Index Index;
    typedef typename TPartitionedPoints::SubDomainIter SubDomainIter;
    typedef typename TPartitionedPoints::Cell Cell;

    typedef typename SubDomainIter::pointer pointer;
    typedef typename SubDomainIter::reference reference;


  public:

    SphereIter() {}

    //If radius < 0 (default), the cutoff of the
    //'owner' PartitionedPoints object will be used
    //instead.

    SphereIter(TPartitionedPoints& owner,
	       TPositionExtractor _posExtractor, 
	       T_num radius=-1) :
      posExtractor(_posExtractor),
      sdIter(owner,false),
      pCell(0)
    { 
      T_num cutoff = owner.getCutoff();
      if( radius < 0 ) {
	radius = cutoff;
	n_cells = 1;
      }
      else if( radius <= cutoff ) {
	n_cells = 1;
      }
      else {
	n_cells = int(std::ceil(radius/cutoff));
      }
      rad2 = radius*radius;
    }

    // Call this at the beginning of iteration loop
    // over the sphere around the 'pointCenter', e.g.
    // as a first element of a 'for(...)' construct.

    void init(const Point& pointCenter) {
      v = pointCenter;
      Index indCenter;
      TPartitionedPoints& owner = sdIter.getOwner();
      owner._getGrid().indAt(pointCenter,indCenter);
      pCell = &(owner._getGrid().getGridArray()(indCenter));
      sdIter.init(indCenter,n_cells);
    }

    // Call as a second element of a 'for(...)' construct

    bool not_end() {
      for( ; sdIter.not_end(); sdIter.next() ) {
	pObject = &(*sdIter);
	dist2 =  blitz_ext::dotSelf(v - posExtractor(*pObject));
	if( dist2 <= rad2 ) 
	  {
	    return true;
	  }
      }
      return false;
    }

    // Alternative call as a second element of a 'for(...)' construct,
    // which does not check the distance but just returns true if
    // the next hashed element is found (that element is made available
    // as the result of dereferencing operator.
    // Works in tandem with 'check_dist()' that actually computes
    // the distance, and that must be called after the call to this
    // method.

    bool not_end_lazy() {
      if( sdIter.not_end() ) {
	pObject = &(*sdIter);
	return true;
      }
      return false;
    }

    // See description of 'not_end_lazy()'

    bool check_dist() {
	dist2 =  blitz_ext::dotSelf(v - posExtractor(*pObject));
	if( dist2 <= rad2 ) 
	  {
	    return true;
	  }
	return false;
    }

    // Call as a third element of a 'for(...)' construct
      
    void next() {
      sdIter.next();
    }


    // Dereferencing this iterator will return a reference 
    // to current PointData object.

    reference operator*() const {
      return *pObject;
    }

    // Pointer dereferencing.

    pointer operator->() const { return pObject; }


    // Power of two for the last computed distance,
    // in other words, the distance^2 between
    // 'pointCenter' provided as an argument to 'init()'
    // and the object returned by 'operator*()'.

    T_num r2() const {
      return dist2;
    }

    // Insert PointData into current neighborhood.
    // Position of this PointData object must be
    // the same as position used in the last call
    // to init() method, otherwise the state
    // of the parent PartitionedPoints object
    // can become incorrect.
    // This way, the call is semantically equivalent
    // to the call to insert() method of the
    // parent PartitionedPoints object, but saves
    // the grid cell lookup. In a typical application
    // of finding neighbors within a set of points,
    // the insertion of the 'pointCenter' is 
    // required right after the loop
    // over the already inserted points. Thus, the cell value
    // calculated within the last init(pointCenter) call can
    // be safely reused.

    void insert(const PointData& x) {
      pCell->push_back(x);
    }

  protected:

    //Object which when called with PointData argument,
    //will return const Point& for the position
    //of that PointData
    TPositionExtractor posExtractor;

    //Sphere radius^2, this radius can be 
    //either larger or smaller than the
    //cutoff distance for the parent PartitionedPoints
    //object
    T_num rad2;

    //Number of PartitionedPoints cells to go through when searching
    //for the objects with a given radius (above).
    //This is computed once in ctor,
    //and then used in each call to init()
    int n_cells;

    //Pointer to the last found PointData object within the sphere.
    //This is the object that this iterator 'points' to.
    pointer pObject;

    //Distance^2 computed for the last found
    //object pObject. It is often required by the calling code,
    //so we cache it here.
    T_num dist2;

    //Point for which this object was created.
    //We need to cache it here because it is used
    //repeatedly by not_end() method
    Point v;

    //Our internal instance of subdomain iterator
    SubDomainIter sdIter;

    //Grid cell found for Point v above;
    Cell *pCell;

  }; // class SphereIter


  //Some of classes which implement PositionExtractor concept
  //to use with SphereIter in typical applications.


#include "PRODDL/Common/common_types.hpp"


  // This helper class will be used for SphereIter
  // when PointData is int index into 1-d Blitz::Array
  // of points in N-d space.

  template<typename T_num,int N_dim>
  struct PositionExtractorByArrayIndex {

    typedef typename PRODDL::common_types::point_type<T_num,N_dim>::Type Point;
    typedef typename PRODDL::common_types::num_vector_of_points_type<T_num,N_dim>::Type Points;

    PositionExtractorByArrayIndex(): pPoints(0) {}

    PositionExtractorByArrayIndex(const Points& _points):
      pPoints(&_points)
    { }

    const Point& operator()(int o) const {
      return (*pPoints)(o);
    }

  protected:
    
    const Points *pPoints;

  };  // struct PositionExtractorByArrayIndex


  // Datatype to use as PointData in PartitionedPoints class.
  // Stores the position and an opaque index into something.


  template<typename T_num,int N_dim>
  struct PointAndIndex {

    typedef typename PRODDL::common_types::point_type<T_num,N_dim>::Type Point;

    Point point;
    int index;

    PointAndIndex(const Point& _point, int _index):
      point(_point),
      index(_index)
    {}

    PointAndIndex() {}

  };  // struct PointAndIndex


  // PositionExtractor implementation for PointAndIndex class

  template<typename T_num,int N_dim>
  struct PositionExtractorFromPointAndIndex {

    typedef typename PRODDL::common_types::point_type<T_num,N_dim>::Type Point;

    const Point& operator()(const PointAndIndex<T_num,N_dim>& o) const {
      return o.point;
    }
  };  // struct PositionExtractorFromPointAndIndex



 } } } // namespace PRODDL { namespace PRODDL { namespace Geom { namespace Points {

#endif // AT_PARTPOINTS_H__
