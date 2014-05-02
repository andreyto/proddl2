//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_GRID_GRID_H__
#define AT_GRID_GRID_H__

//N-dimensional self-extendable rectangular grid that
//can be indexed by spatial coordinates.
//Applications include histograms,
//projections and geometric hashing for
//fast pairwise neighbors lookup.

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Grid/gridcoord.hpp"

#include "PRODDL/Grid/gridexcept.hpp"

#include <blitz/array.h>

#include "PRODDL/Blitz/bzalgor.hpp"

#include "PRODDL/Blitz/bztinyvec.hpp"

#include "PRODDL/Blitz/bzatio.hpp"

#include "PRODDL/Common/logger.hpp"

#include <functional>

namespace PRODDL { namespace Grid {

	//grid object which allowes indexing N-dim array with spatial
	//coord values. at() member
	//function checks for out-of-boundary indexes and  returns
	//reference to dummy variable for out-of-domain indices. at_expand()
	//checks for out-of-bound indices and resizes grid to accomodate
	//them. Grid is expanded by twice the largest size that would include the
	//new index, to make log() number of expansions for the linear sequence of input
	//indices.
	template<int N_dim,class T_val = int,class T_coord = double> class Grid {
	public:
		typedef Geometry::UniformCubicGeometry<N_dim,T_coord> Geom;
		typedef blitz::Array<T_val,N_dim> GridArray;
		typedef blitz::TinyVector<T_coord,N_dim> SpatialCoord;
		typedef blitz::TinyVector<int,N_dim> LogicalCoord;
		typedef Grid<N_dim,T_val,T_coord> _Self;
		typedef blitz::TinyVector<SpatialCoord,2> SpatialDomain;
		typedef blitz::TinyVector<LogicalCoord,2> LogicalDomain;
		typedef T_val T_value;
		enum { N_rank = N_dim };
		//indices of LogicalDomain corresponding to lower and upper bounds:
		enum bounds_ind { lBound = 0, uBound = 1 };
		//flag that defines what to do if index is out-of logical domain
		//domain_nothing - do nothing, expect index must be valid
		//domain_die     - abort
		//domain_expand  - resize the array so that index becomes valid
		//domain_subst   - return the reference to the dummy variable
		enum domain_act_enum { domain_nothing = 0x0, domain_die = 0x2,
			domain_expand = 0x4, domain_subst = 0x8, domain_crop = 0x10 };
		typedef int domain_actions;
		typedef blitz::TinyVector<domain_actions,N_dim> LogicalCoordAct;
		typedef blitz::TinyVector<LogicalCoordAct,2> LogicalDomainAct;
	protected:
		//maps spatial coordinates to array indices
		Geom geom;
		//array with data
		GridArray grid;
		//spatial domain passed to the constructor, can be updated
		//during expansion
		SpatialDomain spa_domain_ini;
		//logical domain i.e. min and max array indices
		LogicalDomain log_domain;
		//reset or not grid with ini_val in ctor and resize():
		bool reset;
		T_value ini_val;
		//dummy variable - reference to it is returned by at() for
		//out-of-domain arguments:
		T_value _x_subst;
		//flag that says what to do if logical domain is violated,
		//hold separate value for lower and upper ends of each dimension.
		//values must be taken from domain_actions enum type
		LogicalDomainAct domain_act;
	public:
		//default constructor:
		Grid() {}

		virtual ~Grid() {}

		//lbound_s and ubound_s - min and max values of spatial coordinates
		//(inclusive); step - spatial coordinate step; reset - if true,
		//initialize grid elements with ini_val.

		Grid(SpatialCoord lbound_s, SpatialCoord ubound_s,
			SpatialCoord step,
			bool reset=true,T_value ini_val=T_value(),
			LogicalCoord lbound = 1) {
				init(lbound_s,ubound_s,step,reset,ini_val,lbound);
		}

		Grid(Geom _geom, LogicalCoord _log_size,
			bool reset=true,T_value ini_val=T_value())
		{
			init(_geom,_log_size,reset,ini_val);
		}

		Grid(SpatialCoord step, LogicalCoord _log_size,
			bool reset=true,T_value ini_val=T_value())
		{
			init(step,_log_size,reset,ini_val);
		}


		// same semantics as above, but data will be referenced
		// from external array 'arr'.
		// The management of data array when grid has to be resized is delegated 
		// to makeNewDataArray() and freeOldDataArray() virtual functions.
		// In default implementation, if array gets resized, the new, independent copy
		// will be created and the connection to the original external array
		// will be lost. Thus, if you use default implementation and want to change data in the external
		// array through this object, do not use 'domain_act' policy
		// that can resize the data (default policy is to raise error on
		// out-of-bound indices in at() function).
		// Note: 'reset' flag will be ignored in ctor (i.e. the initial content
		// of the 'arr' will be always preserved). 'reset' will be used if grid
		// is resized later.

		Grid(SpatialCoord lbound_s, SpatialCoord ubound_s,
			SpatialCoord step,GridArray arr,
			bool reset=true,T_value ini_val=T_value()) {
				init(lbound_s,ubound_s,step,arr,reset,ini_val);
		}

		Grid(Geom _geom, LogicalCoord _log_size,
			GridArray arr,
			bool reset=true,T_value ini_val=T_value())
		{
			init(_geom,_log_size,arr,reset,ini_val);
		}

		Grid(SpatialCoord step, LogicalCoord _log_size,
			GridArray arr,
			bool reset=true,T_value ini_val=T_value())
		{
			init(step,_log_size,arr,reset,ini_val);
		}


		//copy constructor: uses default copy semantics of
		// blitz::Array - data array will reference same memory
		// as source object. If you need to get deep copy, use
		// makeUnique()
		Grid(const Grid& g) : geom(g.geom), grid(g.grid),
			spa_domain_ini(g.spa_domain_ini), log_domain(g.log_domain),
			reset(g.reset),ini_val(g.ini_val), _x_subst(g._x_subst),
			domain_act(g.domain_act) {

		}

		//copy assignment, same semantics as copy ctor above
		_Self& operator= (const _Self& x) {

			if(this != &x) {
				copyButData(x);
				grid.reference(x.grid);
			}
			return *this;
		}


		// assignment of a given T_value to entire grid
		_Self& operator= (const T_value& x) {
			grid = x;
			return *this;
		}

		// make so that this object exclusively owns all it's memory
		// (concerns only grid data array at this time)

		void makeUnique() { 

			GridArray grid_new = makeNewDataArray(LogicalDomain(grid.lbound(),grid.ubound()));
			grid_new = grid;
			grid.reference(grid_new);
			freeOldDataArray();
		}

		//does the actual initialization job - called from ctor and can
		//be called on it's own to rebuild the grid - old content is lost.
		//The ends of range (lbound_s,ubound_s) are excluded.
		void init(SpatialCoord lbound_s, SpatialCoord ubound_s,
			SpatialCoord step,
			bool reset_=true,T_value ini_val_=T_value(),
			LogicalCoord lbound = 1) 
		{

			initButData(lbound_s,ubound_s,step,reset_,ini_val_,lbound);
			//   grid.reference(GridArray(log_domain(lBound),
			//			    LogicalCoord(log_domain(uBound)-log_domain(lBound)+1)));
			grid.reference(makeNewDataArray(log_domain));
			freeOldDataArray();
			if(reset) grid = ini_val;
		}

		// A version that takes as arguments Geom object and size of logical array

		void init(Geom _geom, LogicalCoord _log_size,
			bool reset_=true,T_value ini_val_=T_value())
		{

			initButData(_geom,_log_size,reset_,ini_val_);
			//   grid.reference(GridArray(log_domain(lBound),
			//			    LogicalCoord(log_domain(uBound)-log_domain(lBound)+1)));
			grid.reference(makeNewDataArray(log_domain));
			freeOldDataArray();
			if(reset) grid = ini_val;
		}

		// a version that takes spatial step and  logical array size as arguments,
		// and creates a grid centered around spatial zero point

		void init(SpatialCoord step, LogicalCoord _log_size,
			bool reset_=true,T_value ini_val_=T_value())
		{

			initButData(step,_log_size,reset_,ini_val_);
			grid.reference(makeNewDataArray(log_domain));
			freeOldDataArray();
			if(reset) grid = ini_val;
		}

		// same as init() above, but uses externaly supplied 
		// array for data. In that case,regardless of what 'reset_' is,
		// the array data will not be overwritten with ini_val_ during
		// construction. However, 'reset' flag will be taken into account
		// if grid is resized later.
		void init(SpatialCoord lbound_s, SpatialCoord ubound_s,
			SpatialCoord step, GridArray arr,
			bool reset_=true,T_value ini_val_=T_value())
		{


			initButData(lbound_s,ubound_s,step,reset_,ini_val_,arr.lbound());
			if( ! all(arr.ubound() >= log_domain(uBound)) ) throw grid_sanity_error("External array is too small");
			grid.reference(useExternalDataArray(arr));
			freeOldDataArray();
			//if(reset) grid = ini_val;
		}

		void init(Geom _geom, LogicalCoord _log_size,
			GridArray arr,
			bool reset_=true,T_value ini_val_=T_value())
		{


			initButData(_geom,_log_size,reset_,ini_val_);
			if( ! all(arr.ubound() >= log_domain(uBound)) ) throw grid_sanity_error("External array is too small");
			grid.reference(useExternalDataArray(arr));
			freeOldDataArray();
			//if(reset) grid = ini_val;
		}

		// a version that takes spatial step and  logical array size as arguments,
		// and creates a grid centered around spatial zero point

		void init(SpatialCoord step, LogicalCoord _log_size,
			GridArray arr,
			bool reset_=true,T_value ini_val_=T_value())
		{


			initButData(step,_log_size,reset_,ini_val_);
			if( ! all(arr.ubound() >= log_domain(uBound)) ) throw grid_sanity_error("External array is too small");
			grid.reference(useExternalDataArray(arr));
			freeOldDataArray();
			//if(reset) grid = ini_val;
		}


	protected:

		//Initializes all except data array - must be called from
		//variants of init() function
		//The ends of range (lbound_s,ubound_s) are excluded.
		void initButData(SpatialCoord lbound_s, SpatialCoord ubound_s,
			SpatialCoord step,
			bool reset_=true,T_value ini_val_=T_value(),
			LogicalCoord lbound = 1) {
				spa_domain_ini = lbound_s,ubound_s;
				reset = reset_;
				ini_val = ini_val_;
				LogicalCoord min_lbound = 1;
				//ATALWAYS(all(min_lbound <= lbound),"lbound must be >=1, see \"atgridcoord.h\" for explanation");
				geom = Geom(step,lbound_s,lbound);
				LogicalCoord ubound = geom.toLogical(ubound_s);
				log_domain = 
					lbound, 
					ubound;
				if( ! all(log_domain(lBound) <= log_domain(uBound)) ) throw grid_sanity_error("lBound is more than uBound");
				_x_subst = ini_val;
				//set domain violation action to default value (abort)
				domain_act = domain_die;
		}


		// a version that takes Geom object and logical array size as arguments

		void initButData(Geom _geom, LogicalCoord _log_size, 
			bool reset_=true,T_value ini_val_=T_value()) {
				geom = _geom;
				log_domain = geom.logicalZero(),(_log_size + geom.logicalZero() - 1);
				spa_domain_ini = geom.spatialZero(),geom.toSpatial(log_domain(uBound));
				reset = reset_;
				ini_val = ini_val_;
				_x_subst = ini_val;
				//set domain violation action to default value (abort)
				domain_act = domain_die;
		}

		// a version that takes spatial step and  logical array size as arguments,
		// and creates a grid centered around spatial zero point

		void initButData(SpatialCoord step, LogicalCoord _log_size, 
			bool reset_=true,T_value ini_val_=T_value()) {

				initButData(Geom(step,_log_size),_log_size,reset_,ini_val_);
		}


		void copyButData(const _Self& x) {
			if(this != &x) {
				geom = x.geom;
				spa_domain_ini = x.spa_domain_ini;
				log_domain = x.log_domain;
				reset = x.reset;
				ini_val = x.ini_val;
				_x_subst = x._x_subst;
				domain_act = x.domain_act;
			}
		}

		// The following two virtual functions provide the way to have data array allocated
		// by means other that blitz::Array allocator in derived classes. The data
		// must be, however, in storage format compatible with that of Blitz, so that
		// 'grid' member can be created as referencing preexisting data (neverDeleteData option
		// in Array() ctor). Calls to makeNewDataArray() and freeOldDataArray() must form a sandwitch,
		// and the code in between may have access to both new and old data. See how they are
		// used inside resizeAndPreserve() or makeUnique() methods.

		// Create new data and return GridArray that references it
		// In derived class implementation, it is possible that the returned array 
		// or current 'grid' member 
		// do not own their data, and the old data will be kept alive until
		// 'freeOldDataArray()' is called.
		virtual GridArray makeNewDataArray(const LogicalDomain& ldom_new) {
			return GridArray(ldom_new(lBound),ldom_new(uBound) - ldom_new(lBound) + 1);
		}

		// This method should be called instead of makeNewDataArray when the grid is
		// referencing externally supplied data. See PyGrid for the non-trivial
		// implementation of this virtual method.

		virtual GridArray useExternalDataArray(GridArray arr) {
			return arr;
		}

		// Makes old data to go - must be called after each call to makeNewDataArray().
		// No GridArray managed by this object can reference old data at this time.
		// In particular, 'grid' data member must reference the new array returned by
		// previous call to makeNewDataArray() before this fucntion can be safely called.
		virtual void freeOldDataArray() {
		}

	public:
		void resizeAndPreserve(const LogicalDomain& ldom_new_rel) {
			LogicalDomain ldom_new_abs(log_domain(lBound),log_domain(lBound)+
				(ldom_new_rel(uBound) - ldom_new_rel(lBound)));
			//new array
			GridArray grid_new = makeNewDataArray(ldom_new_abs);
			if( reset ) grid_new = ini_val;
			//intersection of old domain and new relative domain expressed in old domain coords:
			LogicalDomain ldom_o_nr;
			//max of lbound:
			blitz_ext::max_each(log_domain(lBound),ldom_new_rel(lBound),ldom_o_nr(lBound));
			//min of ubound:
			blitz_ext::min_each(log_domain(uBound),ldom_new_rel(uBound),ldom_o_nr(uBound));
			//there is an intersection if for all components (lbound > ubound):
			if(all(ldom_o_nr(lBound) <=  ldom_o_nr(uBound)) ) {
				//intersection expressed in new absolute domain coords:
				LogicalDomain ldom_na_o;
				//BUG!!! there is a "gotcha" in Blitz: if we just write
				//ldom_na_o = ldom_o_nr - ldom_new_rel(lBound) + ldom_new_abs(lBound);
				//then Blitz will interpret it as
				//ldom_na_o(lBound) = ldom_o_nr(lBound) - ldom_new_rel(lBound)(0) + ldom_new_abs(lBound)(0);
				//ldom_na_o(uBound) = ldom_o_nr(uBound) - ldom_new_rel(lBound)(1) + ldom_new_abs(lBound)(1);
				//if N_dim == 2. Hence, we assign each dimension explicitly:
				ldom_na_o(lBound) = ldom_o_nr(lBound) - ldom_new_rel(lBound) + ldom_new_abs(lBound);
				ldom_na_o(uBound) = ldom_o_nr(uBound) - ldom_new_rel(lBound) + ldom_new_abs(lBound);
				//copy data using subdomains:
				grid_new(blitz::RectDomain<N_dim>(ldom_na_o(lBound),ldom_na_o(uBound))) = 
					grid(blitz::RectDomain<N_dim>(ldom_o_nr(lBound),ldom_o_nr(uBound)));
			}
			//set grid array to new value:
			grid.reference(grid_new);
			//do the same for external data that grid might reference
			freeOldDataArray();
			//update other data members:
			SpatialCoord lbound_s = geom.toSpatialLeft(ldom_new_rel(lBound));
			//update coordinate system:
			geom = Geom(geom.spatialStep(), lbound_s, geom.logicalZero());
			//update logical domain:
			log_domain = ldom_new_abs;
			//reconcile spatial domain with logical domain and coordinate system,
			//ARGUABLE: get the center of the bin for the upper bound to avoid
			//extending logical domain to the next bin:
			spa_domain_ini = SpatialDomain(geom.toSpatialLeft(log_domain(lBound)),
				geom.toSpatial    (log_domain(uBound)));
			//_x_subst does not change
			//reset does not change
			//ini_val does not change
			//
			//sanity checks:
			if( !(blitz::all(geom.logicalZero() == ldom_new_abs(lBound))) ) throw grid_sanity_error("Grid::resizeAndPreserve");
			if( !(blitz::all(geom.toLogical(spa_domain_ini(lBound) + geom.spatialStep()/2) == \
				grid.lbound())) ) throw grid_sanity_error("Grid::resizeAndPreserve");
			if( !(blitz_ext::is_each(geom.toLogical(spa_domain_ini(uBound) - geom.spatialStep()/2), \
				LogicalCoord(grid.ubound()+1),std::less<int>())) ) throw grid_sanity_error("Grid::resizeAndPreserve");
		}

		void resizeAndPreserve(const SpatialCoord& s_new_low, const SpatialCoord& s_new_up ) {
			resizeAndPreserve(LogicalDomain(geom.toLogicalUnlimited(s_new_low),geom.toLogicalUnlimited(s_new_up)));
		}

		//return the reference to the cell corresponding to the
		//given spatial coord. No bounds checking is performed
		//unless ATDEBUG_LEVEL >= 4, in which case @see at() is used.
		T_value& operator () (const SpatialCoord& coord_s) {
			LogicalCoord ind(geom.toLogical(coord_s));
#if defined(ATDEBUG) && ATDEBUG_LEVEL >= 4
			ATALWAYS_BEGIN(isInBounds(ind));
			ATDBGERR << "Domain violation: ind=" << ind << " coord_s=" << coord_s << " for grid structure:\n";
			print_short(ATDBGERR);
			ATALWAYS_END;
#endif
			return grid(ind);
		}

		//we have to return value rather than const reference because
		//blitz::Array does it this way
		T_value operator () (const SpatialCoord& coord_s) const {
			LogicalCoord ind(geom.toLogical(coord_s));
#if defined(ATDEBUG) && ATDEBUG_LEVEL >= 4
			ATALWAYS_BEGIN(isInBounds(ind));
			ATDBGERR << "Domain violation: ind=" << ind << " coord_s=" << coord_s << " for grid structure:\n";
			print_short(ATDBGERR);
			ATALWAYS_END;
#endif
			return grid(ind);
		}

		// get the reference to the grid bin corresponding to
		// the spatial coordinate coord_s. If coord_s is out-of-bound,
		// the grid will be expanded to accomodate coord_s*2 if the 
		// corresponding directions have 'domain_act' set to expand,
		// or reference to dummy varyable _x_subst will be returned if
		// 'domain_act' is set to domain_subst and not a single domain_die
		// value, else exception 'grid_domain_error' will be raised

		T_value& at(const SpatialCoord& coord_s) {
			LogicalCoord ind_res;
			if( indAt(coord_s,ind_res) ) {
				return grid(ind_res);
			}
			return _x_subst;
		}


		// does the main job for at() subscript function above
		// takes spatial coord 'coord_s' and sets 'ind_ret' to the
		// corresponding array index and returns true, if valid ind_ret
		// can be found in agreement with current value of domain_act variable,
		// else returns false and leaves ind_ret unchanged.
		// That function will resize the array if 'domain_expand' condition
		// is encountered.
		// Exceptions: raises exception 'grid_domain_error' if 'domain_die'
		// situation is encountered.

		bool indAt(const SpatialCoord& coord_s,LogicalCoord& ind_ret) {
			LogicalCoord ind = geom.toLogical(coord_s);
			const LogicalCoord &_lbound = log_domain(lBound), &_ubound = log_domain(uBound);
			const LogicalCoordAct &_ldomain_act = domain_act(lBound), 
				&_udomain_act = domain_act(uBound);
			domain_actions domain_act_tot = domain_nothing;

			// TODO: optimize

			// if domain_crop, assign the nearest valid value to the index
			for(int j = 0; j < N_dim; j++) {

				if (ind(j) < _lbound(j)) {
					if ( _ldomain_act(j) & domain_crop ) ind(j) = _lbound(j);
				}

				else if (ind(j) > _ubound(j)) {
					if ( _udomain_act(j) & domain_crop ) ind(j) = _ubound(j);
				}
			}

			//accumulate the state of possible domain violations and corresponding
			//actions into
			//one scalar variable domain_act_tot so that it can be tested
			//for the occurence of each particular case listed in domain_actions enum.
			for(int j = 0; j < N_dim; j++) {
				domain_act_tot |=
					(
					( (ind(j) < _lbound(j)) * _ldomain_act(j) ) |
					( (ind(j) > _ubound(j)) * _udomain_act(j) )
					);
			}
			if(domain_act_tot & domain_die) {
				throw grid_domain_error("ATGRID::Grid domain violation");
			}
			else if(domain_act_tot & domain_expand) {
				//find new logical domain expressed relative to the old domain:
				LogicalCoord _lbound_new = _lbound;
				LogicalCoord _ubound_new = _ubound;
				for(int i = 0; i < N_dim; i++) {
					//make new size in each dimension double of the half the maximum distance from
					//any old bound to the index:
					if(ind(i) < _lbound(i)) _lbound_new(i) = ind(i) - (_ubound(i) - ind(i))/2;
					else if(ind(i) > _ubound(i)) _ubound_new(i) = ind(i) + (ind(i) - _lbound(i))/2;
				}
				//expand the grid:
				resizeAndPreserve(LogicalDomain(_lbound_new,_ubound_new));
				//get index using updated coord system:
				ind = geom.toLogical(coord_s);
			}
			else if(domain_act_tot & domain_subst) return false;
			ind_ret = ind;
			return true;
		}

		//return bounds for the spatial domain passed to
		//the ctor. Contrast that to getSpatialDomain() member function.
		//The difference can arise at the upper bound if grid step 
		//is not a round divisor of the spatial range.
		//Invariant: getIniSpatialDomain()(uBound) <= getSpatialDomain()(uBound) &&
		//getIniSpatialDomain()(uBound) + step > getSpatialDomain()(uBound) &&
		//getIniSpatialDomain()(lBound) == getSpatialDomain()(lBound)
		const SpatialDomain& getIniSpatialDomain() const {
			return spa_domain_ini;
		}

		//return allowed domain for spatial coords
		SpatialDomain getSpatialDomain() const {
			return SpatialDomain(geom.toSpatialLeft(log_domain(lBound)),
				geom.toSpatialRight(log_domain(uBound)));
		}

		const LogicalDomain& getLogicalDomain() const {
			return log_domain;
		}

		LogicalCoord getLogicalShape() const {
			return LogicalCoord(log_domain(uBound) - log_domain(lBound) + 1);
		}

		const SpatialCoord& getSpatialStep() const {
			return geom.spatialStep();
		}

		T_coord minSpatialStep() const {

			return blitz::min(getSpatialStep());

		}

		//return true if spatial coordinate is inside valid grid bounds -
		//i.e. if it will be mapped to valid array cell by operator()
		bool isInBounds(const SpatialCoord& coord_s) const {
			return isInBounds(geom.toLogical(coord_s));
		}

		//return true if logical coordinate is inside valid grid bounds -
		//i.e. if it will be mapped to valid array cell by operator()
		bool isInBounds(const LogicalCoord& coord_l) const {
			return blitz_ext::between(coord_l,
				log_domain(lBound),
				log_domain(uBound));
		}

		//bring logical coordinate to valid grid bounds -
		//to the nearest valid cell if initial value is outside of the grid.
		// Return: true if the coordinate was changed ("cropped")
		bool cropLogicalCoord(LogicalCoord& coord_l) const {
			return blitz_ext::crop(coord_l,
				log_domain(lBound),
				log_domain(uBound));
		}

		//bring logical domain to valid grid bounds -
		//lower end of the argument will become no less than the lower end of grid's domain,
		//upper end of the argument will become no less than the upper end of grid's domain
		//Requirement: domain_l[lBound] <= domain_l[uBound]_
		// Return: true if intersection of 'domain_l' with this grid's domain existed at all.
		// Results: If false is returned, the argument is left is invalid state, otherwise it contains
		// the intersection.
		bool cropSubDomain(LogicalDomain& sub_domain_l) const {
			//max of lbound:
			blitz_ext::max_each(log_domain(lBound),sub_domain_l(lBound),sub_domain_l(lBound));
			//min of ubound:
			blitz_ext::min_each(log_domain(uBound),sub_domain_l(uBound),sub_domain_l(uBound));
			//there is an intersection if for all components (lbound > ubound):
			return all(sub_domain_l(lBound) <=  sub_domain_l(uBound));
		}

		const T_value& getIniVal() const {
			return ini_val;
		}
		const LogicalDomainAct& getDomainAction() const {
			return domain_act;
		}
		void setDomainAction(const LogicalDomainAct& new_dom_act_) {
			domain_act = new_dom_act_;
		}
		void setDomainAction(domain_actions new_dom_act_) {
			domain_act = new_dom_act_;
		}
		void setDomainAction(domain_actions new_dom_act_,int iBound,int iDimension) {
			domain_act(iBound)(iDimension) = new_dom_act_;
		}
		bool getReset() const {
			return reset;
		}
		const T_value& getSubst() const {
			return _x_subst;
		} 
		const Geom& getGeometry() const {
			return geom;
		}
		GridArray& getGridArray() {
			return grid;
		}
		const GridArray& getGridArray() const {
			return grid;
		}
		//output for debugging:
		std::ostream& print(std::ostream& out) {
			print_short(out);
			ATSOUTVAR(grid,out); ATSOUTENDL(out);
			return out;
		}
		// print for debugging skipping 
		// grid data
		std::ostream& print_short(std::ostream& out) {
			ATSOUTVAR(geom.spatialStep(),out); ATSOUTENDL(out);
			ATSOUTVAR(spa_domain_ini,out); ATSOUTENDL(out);
			ATSOUTVAR(log_domain,out); ATSOUTENDL(out);
			ATSOUTVAR(reset,out); ATSOUTENDL(out);
			ATSOUTVAR(ini_val,out); ATSOUTENDL(out);
			ATSOUTVAR(_x_subst,out); ATSOUTENDL(out);
			ATSOUTVAR(domain_act,out); ATSOUTENDL(out);
			return out;
		}
		//output for persistence:
		std::ostream& write(std::ostream& out,bool bAsMatrix = true) {
			//output data members nessesary to
			//restore the object later
			out << "# Spatial Domain: [lower bound point, upper bound point]\n";
			out << spa_domain_ini << std::endl;
			out << "# Spatial Step in each dimension\n";
			out << (geom.spatialStep()) << std::endl;
			//output logical_domain: ubound part of it is redundant information and will
			//be used on reading to make sanity check
			out << "# Logical domain [lower bound index, upper bound index]\n";
			out << log_domain << std::endl;
			out << "# Data (C-array order)\n";
			if(bAsMatrix)
				blitz_ext::PrintAsMatrix(out,grid);
			else
				out << grid;
			out << std::endl;
			out << "# Flag: reset array in ctor\n";
			out << reset << std::endl;
			out << "# Initial value to reset array\n";
			out << ini_val << std::endl;
			out << "# Value of a bin used for out-of-bound indices\n";
			out << _x_subst << std::endl;
			out << "# Flags to define action for out-of-bound indices for lower and upper bound in each dimension\n";
			out << domain_act << std::endl;

			return out;
		}

	protected:

		inline std::istream& read_comment(std::istream& in, std::string& scomment) {
			std::string line_end;
			//discard all till the '#' symbol
			std::getline(in,line_end,'#');
			//the rest of line after '#' is comment
			std::getline(in,scomment);
			return in;
		}

	public:

		//input of the object saved by the call to write():
		std::istream& read(std::istream& in) {
			std::string scomment;

			//expect that each data member is preceeded by comment string, which is skipped
			read_comment(in,scomment);
			//DEBUG:
			//OUTVAR(scomment); OUTENDL;

			SpatialDomain spa_domain_ini_;
			in >> spa_domain_ini_;
			//DEBUG:
			//OUTVAR(spa_domain_ini_); OUTENDL;

			read_comment(in,scomment);
			//DEBUG:
			//OUTVAR(scomment); OUTENDL;

			SpatialCoord step_;
			in >> step_;
			//DEBUG:
			//OUTVAR(step_); OUTENDL;

			read_comment(in,scomment);
			//DEBUG:
			//OUTVAR(scomment); OUTENDL;

			//read the logical domain, use the lbound part of it to initialize
			//this grid object, use ubound part later for sanity check

			LogicalDomain log_domain_;
			in >> log_domain_;

			//initialize itself with the new parameters, do not set the grid content to ini_val_:
			init(spa_domain_ini_(lBound),spa_domain_ini_(uBound),
				step_,false,
				T_value(),log_domain_(lBound));

			read_comment(in,scomment);

			//read the grid content
			in >> grid;

			read_comment(in,scomment);

			//read the rest of data members:
			in >> reset;

			read_comment(in,scomment);

			in >> ini_val;

			read_comment(in,scomment);

			in >> _x_subst;

			read_comment(in,scomment);

			in >> domain_act;

			ATASSERT_BEGIN(blitz_ext::f_equal_to(log_domain,log_domain_))
				ATOUTVAR(log_domain); ATOUTENDL();
			ATOUTVAR(log_domain_); ATOUTENDL();
			ATOUTVAR(LogicalDomain(grid.lbound(),grid.ubound())); ATOUTENDL();
			ATOUTLINE("Current content of this grid is (except the data array):\n");
			print_short(ATDBGOUT);
			ATASSERT_END

				return in;
		}

	};


	// return new _Grid object with it's own copy of subdomain in 'g' defined by
	// low and upper corners 's_low' and 's_up'. The copy is made using resizeAndPreserve()
	// member function of _Grid. Thus, subdomain may be beyond the boundaries of original
	// domain. In that case, the cells without corresponding cells in 'g' will be initialized
	// according to settings of 'reset' and 'ini_val' in 'g'.
	// This function is templatized over general _Grid parameter in order to use
	// it with any derived class of Grid<...>.

	template<class _Grid> _Grid getRegionCopyS(const _Grid& g, 
		const typename _Grid::SpatialCoord& s_low,
		const typename _Grid::SpatialCoord& s_up) {
			_Grid g_sub(g); // references data array of g, the rest is separate copy
			g_sub.resizeAndPreserve(s_low,s_up); // now g_sub has it's own data array of new size
			return g_sub;
	}

} } //namespace PRODDL::Grid

#endif //AT_GRID_GRID_H__
