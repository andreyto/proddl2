//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PAIRDIST_H__
#define AT_PAIRDIST_H__

// Build lists of neighbours for points from two sets

#include "PRODDL/Geom/partpoints.hpp"

#include <blitz/array.h>
#include <blitz/tinyvec2.h>

#include <vector>
#include <unordered_map>

namespace PRODDL { namespace Geom { namespace Points {

  // Given two sets of points, for each set find indices of points
  // within cutoff distance from any point in the other set.
  // Precondition: vv1 and vv2 are two sets of points in space,
  // cutoff is cutoof distance, contents of ind1 and ind2 does not matter.
  // Postcondition: ind1 contains unique sorted indices of elements from vv1,
  // ind2 - same for vv2. Previous content of ind1 and ind2 is erased.

  template<class _T_num, int _n_dim> class Neighbors {
  public:
    typedef _T_num T_num;
    enum { n_dim = _n_dim };
    typedef blitz::TinyVector<T_num,n_dim> vect;
    typedef blitz::Array<vect,1> Vvect;
    typedef blitz::TinyVector<int,2> IPair;
    typedef std::vector<IPair> VIPair;
    typedef std::vector<T_num> fvect;
    typedef std::vector<int> ivect;
    typedef std::vector<ivect> ivect2;
    typedef blitz::Array<int,1> iAvect;
    typedef std::unordered_map<int, int> HashInt;
    typedef std::vector<HashInt> VHashInt;


    //class to describe distance from some point to another
    class PointNeighbor
    {
    public:
      int i; //index of another point in some array
      T_num r; //distance to the point i
      PointNeighbor() {}
      PointNeighbor(int _i,T_num _r): i(_i), r(_r) {}
    };

    typedef std::vector<PointNeighbor> PointNeighbors;
    typedef std::vector<PointNeighbors> ListPointNeighbors;


    // Data object to use as PointData parameter for PartitionedPoints2 class

    typedef PointAndIndex<T_num,n_dim> PointAndIndexT;
    typedef PositionExtractorFromPointAndIndex<T_num,n_dim> PointAndIndexPositionExtractor;


/* Subclassing of PartitionedPoints<> class from partpoints.hpp.
   The purpose is to provide simplified interface and methods which
   are more efficient to call from Python (those that deal with entire
   arrays of points at once).
   The public methods from the base class will also be available,
   for cases where 
   Use (point,index) tuples as PointData template parameter;  
   Add methods to process arrays of input points in one call by returning
   lists of points within cutoff distance from each other.
   Example:
   // Compute "electrostatic" force between two sets of points with cutoff:
   PartitionedPoints2 partPoints(lowerBound,upperBound,cutoff);
   partPoints.insert(pointSet1);
   while(...) {
      partPoints.search(pointSet2,pairIndexes,distancesPower2);
      for(int i =0; i < pairIndexes.size(); i++) {        
           force_i = ( pointSet1(pairIndexes(i)(0)) - pointSet2(pairIndexes(i)(1)) ) / 
	         ( sqrt(distancesPower2(i)) * distancesPower2(i) );
           ....
      }
      move_somehow(pointSet2);
   }
*/

    class PartitionedPoints2 : public PartitionedPoints<T_num,PointAndIndexT,n_dim> {


    protected:

      typedef PartitionedPoints<T_num,PointAndIndexT,n_dim> Base;

    public:

      typedef SphereIter<Base,PointAndIndexPositionExtractor> SphereIter2;

    protected:

      // the first set of points has been inserted already if this variable is set to true
      bool firstSetInserted;

    public:

      // Default ctor creates empty object, which must be initialized later by
      // a call to init(...) method.

      PartitionedPoints2() : firstSetInserted(false) {}

      PartitionedPoints2(const typename Base::Point& lBoundS, 
			 const typename Base::Point& uBoundS, 
			 T_num cutoff) :
	Base(lBoundS,uBoundS,cutoff), firstSetInserted(false) 
	{}

      void init(const typename Base::Point& lBoundS, 
		const typename Base::Point& uBoundS, 
		T_num cutoff) {
	Base::init(lBoundS,uBoundS,cutoff);
	firstSetInserted = false;      
      }

      // Insert initial set of points. Insert can be done only once after the object is
      // intitialized or zapPoints() is called.

      void insert(const Vvect& points) {

	ATALWAYS( ! firstSetInserted, \
		  "PartitionedPoints2::insert() method can be called only once");

	for( int i_point =0; i_point < points.size(); i_point++ ) {
	  const vect& point = points(i_point);
	  //old implementation
	  //getCellIndex(point).insert(PointAndIndexT(point,i_point));
	  //new implementation
	  Base::insert(point,PointAndIndexT(point,i_point));
	}

	firstSetInserted = true;

      }

      // Same as the first version of insert(), but also will return pairs
      // of close points within the inserted set.
      // Use this function to search for close pairs of points within one set.
      // Arguments:
      // Input:
      // Vvect& points - array of points
      // Output:
      // VIPair& indexPairs - array of pairs (index of point from 'points',
      // index of point from previously inserted set)
      // fvect& distanceP2 - array of distance^2 between corresponding pairs
      // from 'indexPairs'
      // Output arrays describe only pairs of points with distance between them
      // below cutoff value. 

      void insert(const Vvect& points,VIPair& indexPairs, fvect& distanceP2) {
	do_search(points,indexPairs,distanceP2,true);
      }


      // Search for close pairs of points between the supplied set of points and already
      // inserted points. New points are not inserted.

      void search(const Vvect& points,VIPair& indexPairs, fvect& distanceP2) {
	do_search(points,indexPairs,distanceP2,false);
      }
      

      void zapPoints() {
	Base::zapPoints();
	firstSetInserted = false;
      }
    

    protected:

      void do_search(const Vvect& points,VIPair& indexPairs, 
		     fvect& distanceP2, 
		     bool do_insert) {

	ATALWAYS( ! ( do_insert && firstSetInserted ), \
		  "PartitionedPoints2::insert() method can be called only once");

	indexPairs.clear();
	distanceP2.clear();

	T_num cutoff2 = Base::getCutoffP2();

	for( int i_point =0; i_point < points.size(); i_point++ ) {

	  const vect& v_i = points(i_point);
	  typename Base::CellIndex cellInd = getCellIndex(v_i);
	  for(typename Base::SubDomainIter iterNeighb = cellInd.getNeighbors(); 
	      iterNeighb.not_end(); 
	      iterNeighb.next()) {
	    const PointAndIndexT& point_ind_j = *iterNeighb;
	    T_num r2 =  blitz_ext::dotSelf(v_i-point_ind_j.point);
	    if( r2 <= cutoff2 ) 
	      { 
		indexPairs.push_back(IPair(i_point,point_ind_j.index));
		distanceP2.push_back(r2);
	      }
	  }

	  if( do_insert ) cellInd.insert(PointAndIndexT(v_i,i_point));

	}

	firstSetInserted = true;

      }


      // implementation using SphereIter
      void do_search_new_impl(const Vvect& points,VIPair& indexPairs, 
		     fvect& distanceP2, 
		     bool do_insert) {

	ATALWAYS( ! ( do_insert && firstSetInserted ), \
		  "PartitionedPoints2::insert() method can be called only once");

	indexPairs.clear();
	distanceP2.clear();

	SphereIter2 iter(*this,PointAndIndexPositionExtractor());

	for( int i_point =0; i_point < points.size(); i_point++ ) {

	  const vect& v_i = points(i_point);

	  for(iter.init(v_i); iter.not_end(); iter.next()) {
    
	    // point_ind_j is always within 'radius' distance from 'pointCenter'
	    const PointAndIndexT& point_ind_j = *iter;
	    T_num r2 = iter.r2();
	    indexPairs.push_back(IPair(i_point,point_ind_j.index));
	    distanceP2.push_back(r2);
    
	  }


	  if( do_insert ) iter.insert(PointAndIndexT(v_i,i_point));

	}

	firstSetInserted = true;

      }

    };


    // Take two sets of points vv1,vv2,
    // return in ind1,ind2 indices of points in each set
    // that have at least one point in other set closer than 'cutoff'.

    static void contactPoints(const Vvect& vv1,
			      const Vvect& vv2,
			      T_num cutoff,
			      ivect& ind1,
			      ivect& ind2);


    //given list with two arrays of points vv1,vv2,
    //build lists of 'contact' points for each point (in ind1,ind2),
    //with upper distance limit of 'cutoff'. For the format of output
    //arrays, see structure of 'PointNeighbors' class.
    //Performance: implementation uses PartitionedPoints object, see
    //its description for performance estimates. In case where points are
    //atoms in globular proteins, computation time is O(N_atoms*cutoff^3)

    static void contactPairs(const Vvect& vv1,
			     const Vvect& vv2,
			     T_num cutoff,
			     ListPointNeighbors& ind1,
			     ListPointNeighbors& ind2);

    //version of 'contactPairs' function where contact pairs of points 'ind' are
    //selected from the same set 'vv'.
    static void contactPairs(const Vvect& vv,
			     T_num cutoff,
			     ListPointNeighbors& ind);

    // Function designed primarily to find the lists of contact residues
    // by atom-atom contact criterion: two residues are in
    // contact if they have at least a specified number of atom-atom contacts
    // Arguments: vv1[N1],vv2[N2] - two sets of points ("atoms")
    // ind_group1[N1],ind_group2[N2] - for each point, the index of this point's group (like 'residue number' field in PDB file)
    // 'cutoff' - point-point distance cutoff value, 'minContact' - minimum number of atom-atom contacts between two groups.
    // Results:
    // igroup_cont1[max_group_index_1+1],igroup_cont2[max_group_index_2+1] - for each group, contains a list of groups from other
    // set which are 'in contact' with this group. Groups are in contact if they have at least 'minContacts' contacts between
    // individual elements. Initial content of 'igroup_cont1' and 'igroup_cont2' is discarded
    // through reallocation.

    static void contactGroupPairs(const Vvect& vv1,
				  const Vvect& vv2,
				  const iAvect& ind_group1,
				  const iAvect& ind_group2,
				  T_num cutoff,
				  int minContacts,
				  ivect2& igroup_cont1,
				  ivect2& igroup_cont2);

    // Same as above, but finds contacts within one set of points 'vv'. Contacts within the same group are ignored.

    static void contactGroupPairs(const Vvect& vv,
				  const iAvect& ind_group,
				  T_num cutoff,
				  int minContacts,
				  ivect2& igroup_cont);

    // This function calls contactGroupPairs() and then aggregates its results into two lists: each list contains indices
    // of groups in one set that have at least one contact group in the other set.
    // Arguments: vv1[N1],vv2[N2] - two sets of points ("atoms")
    // ind_group1[N1],ind_group2[N2] - for each point, the index of this point's group (like 'residue number' field in PDB file)
    // 'cutoff' - point-point distance cutoff value
    // Results:
    // igroup_cont1[],igroup_cont2[]

    static void contactGroups(const Vvect& vv1,
			      const Vvect& vv2,
			      const iAvect& ind_group1,
			      const iAvect& ind_group2,
			      T_num cutoff,
			      int minContacts,
			      ivect& igroup_cont1,
			      ivect& igroup_cont2);


    // Accept set of lists of neighbors for a set of points (as range [first,last) of iterators
    // for type 'PointNeighbors'. Select into output range 'result' same set of lists, only
    // references to neighbors with 'r > r_cut' will be excluded from each list. This is equal to
    // calling 'contactPairs' with this value of 'cutoff', but much faster if lists are already short.

    template<class InpIter, class OutIter> static
    OutIter selectFromContactPairs(InpIter first,InpIter last,OutIter result,T_num cutoff)
    {
      for( ; first != last; ++first, ++result)
	{
	  PointNeighbors new_list;
	  for( typename PointNeighbors::const_iterator p = first->begin(); p != first->end(); ++p )
	    if( p->r <= cutoff ) { new_list.push_back(*p); }
	  *result = new_list;
	}
      return result;
    }


    // Straightforward double loop implementation of 'contactPairs' used for testing other implementations.
    static void contactPairsDoubleLoop(const Vvect& vv1,
				       const Vvect& vv2,
				       T_num cutoff,
				       ListPointNeighbors& ind1,
				       ListPointNeighbors& ind2);

    static void contactPairsDoubleLoop(const Vvect& vv,
				       T_num cutoff,
				       ListPointNeighbors& ind);



    // The implementation that uses SphereIter iterator
    static void contactPairsSphereIter(const Vvect& vv1,
				       const Vvect& vv2,
				       T_num cutoff,
				       ListPointNeighbors& ind1,
				       ListPointNeighbors& ind2);


  }; // class Neighbors

} } } // namespace PRODDL::Geom::Points

#include "PRODDL/Geom/pairdist.cc"

#endif // AT_PAIRDIST_H__
