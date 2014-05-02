//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
// Test geom/pairdist.hpp by comparing results of straightforward double-loop algorithm
// with the actual implementation.

#include "PRODDL/Geom/pairdist.hpp"

#include <vector>

#include <algorithm>

#include "PRODDL/Common/math.hpp"

#include "PRODDL/Common/debug.hpp"

#include <boost/timer.hpp>

#include "PRODDL/Common/bz_vect_ext.hpp"

#include <random/uniform.h>

#include "PRODDL/Testing/exception.hpp"

#include "gtest/gtest.h"

  typedef double T_num;

  enum { n_dim = 3 };

  // values smaller than this are considered negligible in comparisons
  const T_num small_val = 1e-8;

  typedef  PRODDL::Geom::Points::Neighbors<T_num,n_dim> Nbs;

  void generateRandomPoints(const Nbs::vect& center, T_num radius, Nbs::Vvect& points) {
    T_num radius2 = radius*radius;
    ranlib::UniformOpenClosed<T_num> rndgen;
    for(int i = 0; i != points.size(); i++) {
      while(true) {
	Nbs::vect& point = points(i);
	// generate randomly oriented vector within cube [-radius,radius]
	for(int j=0; j < n_dim; j++) point(j) = (rndgen.random() - 0.5)*2*radius;
	// if point is within a sphere of required radius, break out of the loop
	// and process the next point
	T_num r2 = blitz_ext::dotSelf(point);
	if(r2 <= radius2) {
	  point += center;
	  break;
	}
      }
    }
  }

  // compare two ListPointNeighbors objects for equality,
  // throw debug exception if they are not equal. Service
  // function for 'test_contactPairs()'

  struct cmpPN {
    bool operator()(const Nbs::PointNeighbor& x, const Nbs::PointNeighbor& y) {
      return x.i < y.i;
    }
  };

  void cmpListPointNeighbors(Nbs::ListPointNeighbors& ind1,
			     Nbs::ListPointNeighbors& ind2) {

    ATALWAYS(ind1.size() == ind2.size(),"Control array of neighbor lists has different length.");
    for(int i = 0; i != ind1.size(); i++) {
      Nbs::PointNeighbors& nb1_i = ind1[i];
      Nbs::PointNeighbors& nb2_i = ind2[i];
      std::sort(nb1_i.begin(),nb1_i.end(),cmpPN());
      std::sort(nb2_i.begin(),nb2_i.end(),cmpPN());
      ATALWAYS_BEGIN(nb1_i.size() == nb2_i.size())
      {
	ATDBGOUT << "Control neighbor list for the point has different length.\n";
	ATOUTVAR(i); ATOUTVAR(nb1_i.size()); ATOUTVAR(nb2_i.size()); ATOUTENDL();
      }
      ATALWAYS_END
      for(int j = 0; j != nb1_i.size(); j++) {
	ATALWAYS_BEGIN(nb1_i[j].i == nb2_i[j].i) 
	  {
	    ATDBGOUT << "Control index of neighbor point is different:\n";
	    ATOUTVAR(j); ATOUTVAR(nb1_i[j].i); ATOUTVAR(nb2_i[j].i); ATOUTENDL();
	    for(int j_n = 0; j_n != nb1_i.size(); j_n++) {
	      ATOUTVAR(nb1_i[j_n].i); ATOUTVAR(nb2_i[j_n].i); ATOUTENDL();
	    }
	  }
	ATALWAYS_END
	//T_num r1 = nb1_i[j].r;
	//T_num r2 = nb2_i[j].r;
	ATALWAYS(PRODDL::Math::are_equal_delta(nb1_i[j].r,nb2_i[j].r,small_val),"Control distance to neighbor point is different.");
      }

    }
  }

TEST(PairdistTest, ContactPairs) {

    const int n_points1 = 1000, n_points2 = 1000;

    // setup parameters for random distributions of points -
    // two partly overlapping spheres
    const T_num cutoff = 1.;
    const T_num distrRadius = cutoff*5;
    const T_num center1 = - distrRadius*1/3;
    const T_num center2 = - center1;
    Nbs::vect vCenter1(center1,0,0), vCenter2(center2,0,0);

    Nbs::Vvect vv1(n_points1), vv2(n_points2);

    // initialize points from random number generator
    generateRandomPoints(vCenter1,distrRadius,vv1);
    generateRandomPoints(vCenter2,distrRadius,vv2);

    typedef std::vector<Nbs::ListPointNeighbors> Vlpn;

    // create two sets for output parameters to call 'contactPairs()': 
    // indMain - for actual implementation, indControl - for control call to
    // contactPairsDoubleLoop()

    Vlpn indMain(2), indControl(2);

    ATDBGOUT << "Testing Default Algorithm vs Double Loop\n";

    // call both main and control algorithms

    boost::timer t;

    Nbs::contactPairs(vv1,vv2,cutoff,indMain[0],indMain[1]);

    double sec = t.elapsed();

    ATDBGOUT << "Default Algorithm time: " << sec << "\n";

    t.restart();

    Nbs::contactPairsDoubleLoop(vv1,vv2,cutoff,indControl[0],indControl[1]);

    sec = t.elapsed();

    ATDBGOUT << "Double Loop Algorithm time: " << sec << "\n";

    // compair results for main algorithm and control algorithm
    // if the results are not equal, exception will be thrown
    cmpListPointNeighbors(indMain[0],indControl[0]);
    cmpListPointNeighbors(indMain[1],indControl[1]);

    ATDBGOUT << "ok\n";

    ATDBGOUT << "Testing SphereIter vs Double Loop\n";

    // same for implementation that uses SphereIter iterator

    t.restart();

    Nbs::contactPairsSphereIter(vv1,vv2,cutoff,indMain[0],indMain[1]);

    sec = t.elapsed();

    ATDBGOUT << "Sphere Iter Algorithm time: " << sec << "\n";

    // compair results for main algorithm and control algorithm
    // if the results are not equal, exception will be thrown
    cmpListPointNeighbors(indMain[0],indControl[0]);
    cmpListPointNeighbors(indMain[1],indControl[1]);

    ATDBGOUT << "ok\n";

  }

