//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_GEOM_CLUSTER_H__
#define PRODDL_GEOM_CLUSTER_H__

#include "PRODDL/types.hpp"
#include "PRODDL/Geom/partpoints.hpp"

#include <functional>
#include <algorithm>

// Classes to cluster points in space

namespace PRODDL {

// Cluster matrix rows;
// The distance is RMSD between rows.
// Each set is supplied a weight (e.g. negated energy function)
// The algorithm calculates the weighted density D_i of rows
// around each set i within the cutoff distance R.
// User supplies the weight array.
// The algorithm then selects as the representative (a "leader"
// for the first cluster) the row with
// the largets D, and assigns all rows within R from it
// to the first cluster.
// This procedure is repeated for the rows not yet assigned
// to any cluster until there are no more unassigned rows.
//
// The implementation hashes first M_hash columns using a
// grid hash to speed up the neighbor search. The M_hash
// template parameter should be chosen carefully because
// the memory for the grid will grow as (Set_Radius/cutoff)^M_hash.
// The default value of M_hash = 5 was selected as a reasonable
// one for protein docking predictions, based on (100 Angstr box size / 
// 5 Angstr cutoff)^5 => 3200000 hash grid cells, which at 
// ~28 bytes overhead for each hash cell element (size of vector<int>
// with one element minus sizeof(int), see src/test/std_container_size.cpp)
// gives ~ 85MB hash grid size.
// This M_hash can be totally off the bat for other applications.
//
// TODO: well, not really linear below, take time to make more exact
// statement.
// Complexity: linear time, linear space on reasonably distributed
// data (e.g. triplets of unit vectors for proteins in close contact)
// Worst case: quadratic time, linear space when first M_hash 
// coordinates of all rows are within rmsdCutoff from each other. 
// 


template<typename T_num,int N_dim>
  struct PositionExtractorClusterMatrix {

    typedef typename PRODDL::common_types::point_type<T_num,N_dim>::Type Point;
    typedef typename PRODDL::common_types::num_multiarray_type<T_num,2>::Type Matrix;

    PositionExtractorClusterMatrix() {}

    explicit
    PositionExtractorClusterMatrix(const Matrix& _m):
      m(_m)
    { }

    PositionExtractorClusterMatrix& operator= 
      (const PositionExtractorClusterMatrix& other)
    { 
      if( this != &other )
	m.reference(other.m);
      return *this;
    }

    Point operator()(int o) const {
      Point p;
      for(int i = 0; i < N_dim; i++)
	p(i) = m(o,i);
      return p;
    }

  protected:
    
    Matrix m;

  };  // struct PositionExtractorClusterMatrix


// Model of Strict Weak Ordering to use in std::sort()
// to sort array indexes by array values.

template<class NumVectType>
struct CompareByIndex : public std::binary_function<int,int,bool> {

  CompareByIndex() { }

  CompareByIndex(const NumVectType& _data) :
    data(_data)
  { }

  CompareByIndex& operator= (const CompareByIndex& other) {
    if( this != &other ) {
      data.reference(other.data);
    }
    return *this;
  }

  bool operator() (int i, int j) const {

    return data(i) < data(j);

  }

  NumVectType data;

}; // CompareByIndex

template<typename _T_num, int _M_hash=4>
  class ClusterMatrix {

  public:

    typedef _T_num T_num;
    enum { M_hash = _M_hash };

    typedef typename PRODDL::common_types::num_multiarray_type<T_num,2>::Type Matrix;
    typedef Geom::Points::PartitionedPoints<T_num,int,M_hash> PartPoints;
    typedef typename PartPoints::Point PartPoint;
    typedef PositionExtractorClusterMatrix<T_num,M_hash> PosExtractor;
    typedef Geom::Points::SphereIter<PartPoints,PosExtractor> SpIter;
    typedef typename PRODDL::common_types::num_vector_type<T_num>::Type Floats;
    typedef typename PRODDL::common_types::num_vector_type<int>::Type Ints;

  protected:

    Ints leader;
    Floats dens;
    Ints cnt;
    Ints ind_dens;
    Floats dist;

    int n_leaders;

  public:
    
    ClusterMatrix() : n_leaders(0) {}

    void
    cluster(const Matrix& m, const Floats& weights, T_num rmsdCutoff) {

      using namespace blitz;

      // 

      int M_col = m.cols();

      int M_colRest = M_col - M_hash;

      ATLOG_OUT_5(ATLOGVAR(M_col) << ATLOGVAR(M_hash));

      ATLOG_ASSERT_1(M_col >= M_hash);

      T_num rmsdCutoff2 = rmsdCutoff*rmsdCutoff;

      T_num cutoff2 = rmsdCutoff2 * M_col;

      T_num cutoff = std::sqrt(cutoff2);

      // number of sets
      int N_rows = m.rows();

      ATLOG_ASSERT_1(weights.rows() == N_rows);

      // Create a view of the input matrix for the first
      // M_hash columns

      Matrix m_hash(m(Range::all(),Range(0,M_hash-1)));

      // And a view for the rest of the matrix

      Range rColRest(M_hash,toEnd);

      // just some valid Range value if there are no columns left
      // for m_rest. m_rest will not be used anyway

      if(M_colRest == 0) {
	rColRest = Range::all();
      }

      Matrix m_rest(m(Range::all(),rColRest));

      // Density array that will be computed here

      dens.resize(N_rows);

      // Initially, density includes the element itself

      dens = weights;

      cnt.resize(N_rows);

      // Initially, count is 1, meaning element is in its own neigborhood
      
      cnt = 1;

      // Find the bounds, initialize the PartitionedPoints structure
      // and corresponding sphere iterator

      firstIndex i_ind;
      secondIndex j_ind;

      Floats a_min(min(m_hash(j_ind,i_ind),j_ind));
      Floats a_max(max(m_hash(j_ind,i_ind),j_ind));

      // convert bounds to the type required by PartPoints
      // ctor

      PartPoint b_min, b_max;

      for( int i = 0; i < M_hash; i++ ) {
	b_min(i) = a_min(i);
	b_max(i) = a_max(i);
      }

      ATLOG_OUT_5(ATLOGVAR(b_min) << ATLOGVAR(b_max));

      // We create partPoints on the stack, so that
      // the memory will be freed after a call to this
      // method. If the usage pattern requires repeated
      // calls to this method with similarily distributed
      // input, we might have to introduce some policy
      // option into the ctor to keep partPoints between
      // the calls.
      PartPoints partPoints;
      partPoints.init(b_min,b_max,cutoff);
      SpIter sphereIter(partPoints,PosExtractor(m));

      // find the density

      for(int i_row = 0; i_row < N_rows; i_row++) {

	PartPoint pointHash;

	for(int j_col = 0; j_col < M_hash; j_col++) {
	  pointHash(j_col) = m(i_row,j_col);
	}

	for(sphereIter.init(pointHash); sphereIter.not_end(); sphereIter.next()) {
	    int k_row = *sphereIter;

	    T_num r2 = sphereIter.r2();

	    bool inCutoff = true;

	    if( M_colRest > 0 ) {

	      // we already have r^2 for the hashed part of the row,
	      // let's find the same for other columns
	      T_num r2_rest = blitz_ext::dotSelf(m_rest(i_row,Range::all()) - 
						 m_rest(k_row,Range::all()));
	      
	      r2 += r2_rest;

	      // sphereIter.r2() is already < cutoff2,
	      // but we still need to check the total r^2

	      if( r2 >= cutoff2 ) {
		inCutoff = false;
	      }

	    }


	    if( inCutoff ) {

	      dens(i_row) += weights(k_row);
	      dens(k_row) += weights(i_row);

	      cnt(i_row)++;
	      cnt(k_row)++;

	    }
	}
	
	sphereIter.insert(i_row);

      }

      //TMP:
      dens = weights;

      // sort the row index by density

      leader.resize(N_rows);
      ind_dens.resize(N_rows);

      // at first cluster leader and density order
      // indexes point each row to itself

      for(int i_row = 0; i_row < N_rows; i_row++) {
	leader(i_row) = i_row;
	ind_dens(i_row) = i_row;
      }

      // sort density order index
      // blitz-0.7 "stl-compliant" iterators do not work with std::sort
      // (lack of operator- ), so we use raw pointers here.
      // this is safe because ind_dens was allocated locally
      // and therefore is always contigues
      std::sort(ind_dens.dataFirst(),ind_dens.dataFirst()+N_rows,
		CompareByIndex<Floats>(dens));

      // reverse lookup density index - stores
      // in i-th element the place of i-th row
      // in 'ind_dens'
      Ints x_ind_dens(N_rows);

      for(int i_row = 0; i_row < N_rows; i_row++) {
	x_ind_dens(ind_dens(i_row)) = i_row;
      }

      // iterate the rows in decreasing density and 
      // assign cluster leaders

      // 'dist' is the distance from each row to its cluster leader

      dist.resize(N_rows);

      dist = 0;

      n_leaders = 0;

      // we will be going through each pair of neighbors
      // only once, so zap the hash and start inserting
      // points again

      partPoints.zapPoints();

      for(int i_dens = 0; i_dens < N_rows; i_dens++) {

	int i_row = ind_dens(i_dens);

	PartPoint pointHash;

	for(int j_col = 0; j_col < M_hash; j_col++) {
	  pointHash(j_col) = m(i_row,j_col);
	}

	int i_leader = leader(i_row); // must be equal to i_row at this point
	int place_dens = x_ind_dens(i_leader); // top density index of any leader found so far for row
	T_num r2_i = 0; // distance to i_leader

	for(sphereIter.init(pointHash); sphereIter.not_end(); sphereIter.next()) {
	    int k_row = *sphereIter;

	    T_num r2 = sphereIter.r2();

	    bool inCutoff = true;

	    if( M_colRest > 0 ) {

	      T_num r2_rest = blitz_ext::dotSelf(m_rest(i_row,Range::all()) - 
						 m_rest(k_row,Range::all()));

	      r2 += r2_rest;

	      if( r2 >= cutoff2 ) {
		inCutoff = false;
	      }

	    }

	    if( inCutoff ) {

	      // we only inserted cluster leaders into the hash, so
	      // k is a cluster leader already and is within cutoff from i,
	      // and was inserted before i in the order of 
	      // highest density,
	      // thus, k is a candidate for a leader for i.
	      // but we still need to select the one with
	      // the highest density (the lowest ind_dens)
	      // among such candidates
	      // we compare by ind_dens to select the same
	      // leader irrespective of iteration order
	      // in cases when we have equal density for
	      // alternative leaders.

	      int place_dens_k = x_ind_dens(k_row);
	      if( place_dens_k < place_dens ) {
		place_dens = place_dens_k;
		i_leader = k_row;
		r2_i = r2;
	      }

	    }

	}
	
	// if i has not been assigned a leader till this point,
	// its leader index will stay as it was, which is itself
	// thus, i will become a new leader
	// we only insert leaders into the hash
	if( i_row == i_leader ) {
	  sphereIter.insert(i_row);
	  n_leaders++;
	}
	else {
	  leader(i_row) = i_leader;
	  dist(i_row) = std::sqrt(r2_i);
	}
      }
    }


    // Accessor methods for results of the last
    // call to 'cluster()'

    const Ints& getLeader() const {
      return leader;
    }

    const Ints& getCount() const {
      return cnt;
    }

    const Floats& getDensity() const {
      return dens;
    }

    const Ints& getDensityIndex() const {
      return ind_dens;
    }

    const Floats& getDistance() const {
      return dist;
    }

    int numClusters() const {
      return n_leaders;
    }

    // Take output arrays with size equal or greater than
    // the number of clusters returned by numClusters() 
    // and store one record for each cluster.
    // Arguments:
    // Output: 
    // 'rowInd' - index of each cluster leader.
    // 'clustDens' - density for each leader pointed by 'rowInd'
    // Results are sorted from larger to smaller density.
    // Precondition: call to 'cluster()'

    void getClusters(Ints& rowInd, Floats& clustDens) const {
      ATLOG_ASSERT_1(rowInd.rows() >= n_leaders && \
		     clustDens.rows() >= n_leaders);
      int N_rows = ind_dens.rows();
      // go through original rows from high density to low density
      int i_clust = 0;
      for( int i_dens = 0; i_dens < N_rows; i_dens++ ) {
	int i_row = ind_dens(i_dens);
	if( leader(i_row) == i_row ) { // this is a cluster leader
	  ATLOG_ASSERT_1(i_clust < n_leaders);
	  rowInd(i_clust) = i_row;
	  //TMP:
	  clustDens(i_clust) = T_num(cnt(i_row));
	  //clustDens(i_clust) = dens(i_row);
	  i_clust++;
	}
      }
    }

  };

} // namespace PRODDL

#endif // PRODDL_GEOM_CLUSTER_H__
