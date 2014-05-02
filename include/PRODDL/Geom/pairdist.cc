#ifndef  AT_PAIRDIST_CC__
#define  AT_PAIRDIST_CC__

#ifndef AT_PAIRDIST_H__
#  error This file must be included by pairdist.hpp
#endif

#include "PRODDL/Common/bz_vect_ext.hpp"

#include "PRODDL/Blitz/bzarray_iter.hpp"

namespace PRODDL { namespace Geom { namespace Points {

  template<class _T_num, int _n_dim> 
  void Neighbors<_T_num,_n_dim>::contactPoints(const typename Neighbors<_T_num,_n_dim>::Vvect& vv1,
						      const typename Neighbors<_T_num,_n_dim>::Vvect& vv2,
						      _T_num cutoff,
						      typename Neighbors<_T_num,_n_dim>::ivect& ind1,
						      typename Neighbors<_T_num,_n_dim>::ivect& ind2)
  {
    T_num cutoff2 = blitz::pow2(cutoff);
    ind1.clear(); ind2.clear();

    //pairwise distance calculation:
    for(int i = 0; i != vv1.size(); i++) for(int j = 0; j != vv2.size(); j++)
      if( blitz_ext::dotSelf(vv1(i)-vv2(j)) <= cutoff2 ) { ind1.push_back(i); ind2.push_back(j); }

    //make selected indices unique:
    sort(ind1.begin(),ind1.end());
    ind1.erase(unique(ind1.begin(),ind1.end()),ind1.end());
    sort(ind2.begin(),ind2.end());
    ind2.erase(unique(ind2.begin(),ind2.end()),ind2.end());
  }

  // Find pairs of points from two sets vv1 and vv2 such that points from different sets are within cutoff distance from each other.
  // Straightforward double loop implementation for testing other implementations.

  template<class _T_num, int _n_dim> 
  void Neighbors<_T_num,_n_dim>::contactPairsDoubleLoop(
	  const typename Neighbors<_T_num,_n_dim>::Vvect& vv1,
           const typename Neighbors<_T_num,_n_dim>::Vvect& vv2,
           _T_num cutoff,
           typename Neighbors<_T_num,_n_dim>::ListPointNeighbors& ind1,
	   typename Neighbors<_T_num,_n_dim>::ListPointNeighbors& ind2)
  {

    T_num cutoff2 = blitz::pow2(cutoff);
    //reset output array
    ind1.clear();
    ind1.resize(vv1.size());
    ind2.clear();
    ind2.resize(vv2.size());
    for(int i = 0; i != vv1.size(); i++) 
      {
	PointNeighbors& nb_i = ind1[i];
	for(int j = 0; j != vv2.size(); j++)
	  {
	    T_num r2 =  blitz_ext::dotSelf(vv1(i)-vv2(j));
	    if( r2 <= cutoff2 ) 
	      { 
		T_num r = sqrt(r2);
		PointNeighbors& nb_j = ind2[j];
		nb_i.push_back(PointNeighbor(j,r));
		nb_j.push_back(PointNeighbor(i,r));
	      }
	  }
      }
  }


  // Find pairs of points from set vv such that points are within cutoff distance from each other.
  // Straightforward double loop implementation for testing other implementations.

  template<class _T_num, int _n_dim> 
  void Neighbors<_T_num,_n_dim>::contactPairsDoubleLoop(const typename Neighbors<_T_num,_n_dim>::Vvect& vv,
							_T_num cutoff,
							typename Neighbors<_T_num,_n_dim>::ListPointNeighbors& ind)
  {

    T_num cutoff2 = blitz::pow2(cutoff);
    //reset output array
    ind.clear();
    ind.resize(vv.size());
    for(int i = 0; i != vv.size(); i++) 
      {
	PointNeighbors& nb_i = ind[i];
	for(int j = 0; j < i; j++)
	  {
	    T_num r2 =  blitz_ext::dotSelf(vv(i)-vv(j));
	    if( r2 <= cutoff2 ) 
	      { 
		T_num r = sqrt(r2);
		PointNeighbors& nb_j = ind[j];
		nb_i.push_back(PointNeighbor(j,r));
		nb_j.push_back(PointNeighbor(i,r));
	      }
	  }
      }
  }


  template<class _T_num, int _n_dim> 
  void Neighbors<_T_num,_n_dim>::
  contactPairsSphereIter(const typename Neighbors<_T_num,_n_dim>::Vvect& vv1,
			 const typename Neighbors<_T_num,_n_dim>::Vvect& vv2,
			 _T_num cutoff,
			 typename Neighbors<_T_num,_n_dim>::ListPointNeighbors& ind1,
			 typename Neighbors<_T_num,_n_dim>::ListPointNeighbors& ind2)
  {

    typedef PositionExtractorByArrayIndex<_T_num,_n_dim> IndexPositionExtractor;

    //reset output array
    ind1.clear();
    ind1.resize(vv1.size());
    ind2.clear();
    ind2.resize(vv2.size());
    if(vv1.size() == 0 || vv2.size() == 0) return; // so that bracket() below will work guaranteed
    typedef PartitionedPoints<T_num,int,n_dim> PartPoints;
    typedef typename PartPoints::Point Point;
    typedef SphereIter<PartPoints,IndexPositionExtractor> SphereIter;
    Point lbound_s = vv1(0);
    Point ubound_s = lbound_s;
    //TODO: currently grid size is not restricted, 
    //may get huge for distant molecules, for instance.
    for(int i = 0; i < vv1.size(); i++) {
      blitz_ext::bracket(vv1(i),lbound_s,ubound_s);
    }
    for(int i = 0; i < vv2.size(); i++) {
      blitz_ext::bracket(vv2(i),lbound_s,ubound_s);
    }
    PartPoints partPoints(lbound_s,ubound_s,cutoff);
    for(int i = 0; i != vv2.size(); i++) 
      {
	partPoints.insert(vv2(i),i);
      }

    SphereIter iter(partPoints,IndexPositionExtractor(vv2));

    for(int i = 0; i != vv1.size(); i++) 
      {
	PointNeighbors& nb_i = ATTAKE_DBG(ind1,i);
	const vect& v1_i = vv1(i);
	for(iter.init(v1_i); iter.not_end(); iter.next()) {
	  int j = *iter;
	  T_num r2 =  iter.r2();
	  T_num r = sqrt(r2);
	  PointNeighbors& nb_j = ATTAKE_DBG(ind2,j);
	  nb_i.push_back(PointNeighbor(j,r));
	  nb_j.push_back(PointNeighbor(i,r));
	}
      }
  }



  template<class _T_num, int _n_dim> 
  void Neighbors<_T_num,_n_dim>::
  contactPairs(const typename Neighbors<_T_num,_n_dim>::Vvect& vv1,
	       const typename Neighbors<_T_num,_n_dim>::Vvect& vv2,
	       _T_num cutoff,
	       typename Neighbors<_T_num,_n_dim>::ListPointNeighbors& ind1,
	       typename Neighbors<_T_num,_n_dim>::ListPointNeighbors& ind2)
  {

    T_num cutoff2 = blitz::pow2(cutoff);
    //reset output array
    ind1.clear();
    ind1.resize(vv1.size());
    ind2.clear();
    ind2.resize(vv2.size());
    if(vv1.size() == 0 || vv2.size() == 0) return; // so that bracket() below will work garanteed
    typedef PartitionedPoints<T_num,int,n_dim> PartPoints;
    typedef typename PartPoints::Point Point;
    Point lbound_s = vv1(0);
    Point ubound_s = lbound_s;
    //TODO: currently grid size is not restricted, 
    //may get huge for distant molecules, for instance.
    for(int i = 0; i < vv1.size(); i++) {
      blitz_ext::bracket(vv1(i),lbound_s,ubound_s);
    }
    for(int i = 0; i < vv2.size(); i++) {
      blitz_ext::bracket(vv2(i),lbound_s,ubound_s);
    }
    PartPoints partPoints(lbound_s,ubound_s,cutoff);
    for(int i = 0; i != vv2.size(); i++) 
      {
	partPoints.getCellIndex(vv2(i)).insert(i);
      }
    for(int i = 0; i != vv1.size(); i++) 
      {
	PointNeighbors& nb_i = ATTAKE_DBG(ind1,i);
	const vect& v1_i = vv1(i);
	typename PartPoints::CellIndex cellInd = partPoints.getCellIndex(v1_i);
	for(typename PartPoints::SubDomainIter iterNeighb = cellInd.getNeighbors(); 
	    iterNeighb.not_end(); 
	    iterNeighb.next()) {
	  int j = *iterNeighb;
	  T_num r2 =  blitz_ext::dotSelf(v1_i-vv2(j));
	  if( r2 <= cutoff2 ) 
	    { 
	      T_num r = sqrt(r2);
	      PointNeighbors& nb_j = ATTAKE_DBG(ind2,j);
	      nb_i.push_back(PointNeighbor(j,r));
	      nb_j.push_back(PointNeighbor(i,r));
	    }
	}
      }
  }


  template<class _T_num, int _n_dim> 
  void Neighbors<_T_num,_n_dim>::contactPairs(const typename Neighbors<_T_num,_n_dim>::Vvect& vv,
					      _T_num cutoff,
					      typename Neighbors<_T_num,_n_dim>::ListPointNeighbors& ind)
  {

    T_num cutoff2 = blitz::pow2(cutoff);
    //reset output array
    ind.clear();
    ind.resize(vv.size());
    if(vv.size() == 0) return; // so that bracket() below will work guaranteed
    typedef PartitionedPoints<T_num,int,n_dim> PartPoints;
    typedef typename PartPoints::Point Point;
    Point lbound_s = vv(0);
    Point ubound_s = lbound_s;
    //TODO: currently grid size is not restricted, may get huge for distant molecules, for instance.
    for(int i = 0; i < vv.size(); i++) {
      blitz_ext::bracket(vv(i),lbound_s,ubound_s);
    }
    PartPoints partPoints(lbound_s,ubound_s,cutoff);
    for(int i = 0; i != vv.size(); i++) 
      {
	PointNeighbors& nb_i = ind[i];
	const vect& v_i = vv(i);
	typename PartPoints::CellIndex cellInd = partPoints.getCellIndex(v_i);
	for(typename PartPoints::SubDomainIter iterNeighb = cellInd.getNeighbors(); 
	    iterNeighb.not_end(); 
	    iterNeighb.next()) {
	  int j = *iterNeighb;
	  T_num r2 =  blitz_ext::dotSelf(v_i-vv(j));
	  if( r2 <= cutoff2 ) 
	    { 
	      T_num r = sqrt(r2);
	      PointNeighbors& nb_j = ind[j];
	      nb_i.push_back(PointNeighbor(j,r));
	      nb_j.push_back(PointNeighbor(i,r));
	    }
	}
	cellInd.insert(i);
      }
  }



  template<class _T_num, int _n_dim> 
  void Neighbors<_T_num,_n_dim>::contactGroupPairs(const typename Neighbors<_T_num,_n_dim>::Vvect& vv1,
						   const typename Neighbors<_T_num,_n_dim>::Vvect& vv2,
						   const typename Neighbors<_T_num,_n_dim>::iAvect& ind_group1,
						   const typename Neighbors<_T_num,_n_dim>::iAvect& ind_group2,
						   _T_num cutoff,
						   int minContacts,
						   typename Neighbors<_T_num,_n_dim>::ivect2& igroup_cont1,
						   typename Neighbors<_T_num,_n_dim>::ivect2& igroup_cont2)
  {

    ListPointNeighbors ind1, ind2;
    contactPairs(vv1,vv2,cutoff,ind1,ind2);
    // zap output arrays
    igroup_cont1.clear(); igroup_cont2.clear();
    // if any of group arrays empty (means vv1 or vv2 must be emty too),
    // return immediately (so that max_element() below is protected)
    if(ind_group1.size() == 0 || ind_group2.size() == 0)
      return;
    // find max group index for each set and resize output arrays to max+1
    int max_group1 = *std::max_element(ind_group1.begin(),ind_group1.end());
    igroup_cont1.resize(max_group1+1);
    int max_group2 = *std::max_element(ind_group2.begin(),ind_group2.end());
    igroup_cont2.resize(max_group2+1);
    // create arrays of hashes to compute number of contacts in each group
    VHashInt vhash1(igroup_cont1.size()),vhash2(igroup_cont2.size()) ;
    // go through point-point contacts and accumulate group-group contacts,
    // possibly not unique    
    for(int i = 0; i != ind1.size(); i++) {
      PointNeighbors& nb_i = ind1[i];
      int i_group1 = ind_group1(i);
      for(int j = 0; j != nb_i.size(); j++) {
	int i_group2 = ind_group2(nb_i[j].i);
	vhash1[i_group1][i_group2]++;
	vhash2[i_group2][i_group1]++;
	//igroup_cont1[i_group1].push_back(i_group2);
	//igroup_cont2[i_group2].push_back(i_group1);
      }
    }
    // make indices of contact groups unique for each group
    for(int i = 0; i != igroup_cont1.size(); i++) {
      ivect& igroup = igroup_cont1[i];
      const HashInt& hashCont = vhash1[i];
      for(HashInt::const_iterator hashIter = hashCont.begin(); hashIter != hashCont.end(); hashIter++) {
	int i_other_group_ind = hashIter->first;
	int i_other_group_ncont = hashIter->second;
	// exclude groups with contacts less than 'minContacts'
	if(i_other_group_ncont >= minContacts)
	  igroup.push_back(i_other_group_ind);
      }
      std::sort(igroup.begin(),igroup.end());
      //igroup.erase(std::unique(igroup.begin(),igroup.end()),igroup.end());      
    }
    for(int i = 0; i != igroup_cont2.size(); i++) {
      ivect& igroup = igroup_cont2[i];
      const HashInt& hashCont = vhash2[i];
      for(HashInt::const_iterator hashIter = hashCont.begin(); hashIter != hashCont.end(); hashIter++) {
	int i_other_group_ind = hashIter->first;
	int i_other_group_ncont = hashIter->second;
	// exclude groups with contacts less than 'minContacts'
	if(i_other_group_ncont >= minContacts)
	  igroup.push_back(i_other_group_ind);
      }
      std::sort(igroup.begin(),igroup.end());
      //igroup.erase(std::unique(igroup.begin(),igroup.end()),igroup.end());      
    }
  }


  template<class _T_num, int _n_dim> 
  void Neighbors<_T_num,_n_dim>::contactGroupPairs(const typename Neighbors<_T_num,_n_dim>::Vvect& vv,
						   const typename Neighbors<_T_num,_n_dim>::iAvect& ind_group,
						   _T_num cutoff,
						   int minContacts,
						   typename Neighbors<_T_num,_n_dim>::ivect2& igroup_cont)
  {

    ListPointNeighbors ind;
    contactPairs(vv,cutoff,ind);
    // zap output array
    igroup_cont.clear();
    // if group array empty (means vv must be emty too),
    // return immediately (so that max_element() below is protected)
    if(ind_group.size() == 0)
      return;
    // find max group index and resize output array to max+1
    int max_group = *std::max_element(ind_group.begin(),ind_group.end());
    igroup_cont.resize(max_group+1);
    // create array of hashes to compute number of contacts in each group
    VHashInt vhash(igroup_cont.size());
    // go through point-point contacts and accumulate group-group contacts,
    // possibly not unique    
    for(int i = 0; i != ind.size(); i++) {
      PointNeighbors& nb_i = ind[i];
      int i_group1 = ind_group(i);
      for(int j = 0; j != nb_i.size(); j++) {
	int i_group2 = ind_group(nb_i[j].i);
	vhash[i_group1][i_group2]++;
	vhash[i_group2][i_group1]++;
	//igroup_cont[i_group1].push_back(i_group2);
	//igroup_cont[i_group2].push_back(i_group1);
      }
    }
    // make indices of contact groups unique for each group
    for(int i = 0; i != igroup_cont.size(); i++) {
      ivect& igroup = igroup_cont[i];
      const HashInt& hashCont = vhash[i];
      for(HashInt::const_iterator hashIter = hashCont.begin(); hashIter != hashCont.end(); hashIter++) {
	int i_other_group_ind = hashIter->first;
	int i_other_group_ncont = hashIter->second;
	// exclude self contacts and groups with contacts less than 'minContacts'
	if(i_other_group_ind != i && i_other_group_ncont >= minContacts)
	  igroup.push_back(i_other_group_ind);
      }
      std::sort(igroup.begin(),igroup.end());
      //igroup.erase(std::unique(igroup.begin(),igroup.end()),igroup.end());      
    }
  }



  template<class _T_num, int _n_dim> 
  void Neighbors<_T_num,_n_dim>::contactGroups(const typename Neighbors<_T_num,_n_dim>::Vvect& vv1,
					       const typename Neighbors<_T_num,_n_dim>::Vvect& vv2,
					       const typename Neighbors<_T_num,_n_dim>::iAvect& ind_group1,
					       const typename Neighbors<_T_num,_n_dim>::iAvect& ind_group2,
					       _T_num cutoff,
					       int minContacts,
					       typename Neighbors<_T_num,_n_dim>::ivect& igroup_cont1,
					       typename Neighbors<_T_num,_n_dim>::ivect& igroup_cont2)
  {

    //ATDBGOUT << "ATDEBUG_CONTACTS: ";
    //ATOUTVAR(ind_group1.size()); ATOUTVAR(ind_group2.size());
    //ATOUTVAR(vv1.size()); ATOUTVAR(vv2.size()); ATOUTENDL();

    // zap output arrays
    igroup_cont1.clear(); igroup_cont2.clear();

    // if any of group arrays empty (means vv1 or vv2 must be emty too),
    // return immediately (so that max_element() below is protected)
    if(ind_group1.size() == 0 || ind_group2.size() == 0)
      return;

    ListPointNeighbors ind1, ind2;
    contactPairs(vv1,vv2,cutoff,ind1,ind2);

    // find max group index for each set and resize output arrays to max+1
    int max_group1 = *std::max_element(ind_group1.begin(),ind_group1.end());
    igroup_cont1.resize(max_group1+1,0);
    int max_group2 = *std::max_element(ind_group2.begin(),ind_group2.end());
    igroup_cont2.resize(max_group2+1,0);

    // go through point-point contacts and accumulate group contacts
    for(int i = 0; i != ind1.size(); i++) {
      PointNeighbors& nb_i = ATTAKE_DBG(ind1,i);
      int i_group1 = ind_group1(i);
      for(int j = 0; j != nb_i.size(); j++) {
	int i_group2 = ind_group2(ATTAKE_DBG(nb_i,j).i);
	ATTAKE_DBG(igroup_cont1,i_group1)++;
	ATTAKE_DBG(igroup_cont2,i_group2)++;
      }
    }

    // store 'in place' indices of those groups that have more than 'minContacts' total contacts

    int j_group = 0;
    for(int i_group = 0; i_group < igroup_cont1.size(); i_group++) {
      if(ATTAKE_DBG(igroup_cont1,i_group) >= minContacts)
	ATTAKE_DBG(igroup_cont1,j_group++) = i_group;
    }
    igroup_cont1.erase(igroup_cont1.begin()+j_group,igroup_cont1.end());

    j_group = 0;
    for(int i_group = 0; i_group < igroup_cont2.size(); i_group++) {
      if(ATTAKE_DBG(igroup_cont2,i_group) >= minContacts)
	ATTAKE_DBG(igroup_cont2,j_group++) = i_group;
    }
    igroup_cont2.erase(igroup_cont2.begin()+j_group,igroup_cont2.end());
  }


} } } // namespace PRODDL::Geom::Points

#endif // AT_PAIRDIST_CC__
