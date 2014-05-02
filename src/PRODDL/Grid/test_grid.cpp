//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Grid/grid.hpp"

#include "PRODDL/Common/debug.hpp"

#include "gtest/gtest.h"

using namespace PRODDL;

typedef Grid::Grid<3,int,double> Grid3D;
typedef Grid3D::SpatialCoord PointS;
typedef Grid3D::LogicalCoord PointL;

TEST(GridTest, All) {

    Grid3D g(PointS(-50,-50,-50),PointS(50,50,50),PointS(20,20,20),true,1,PointL(-50,-50,-50));

	g.at(PointS(-20,-20,-20)) = 2;

	ATALWAYS(g.at(PointS(-20,-20,-20)) == 2,"Unexpected value in grid after assignment");

	g.at(PointS(-20,-20,-20)) += 30;

	ATALWAYS(g.at(PointS(-20,-20,-20)) == 32,"Unexpected value in grid after increment");

	ASSERT_ANY_THROW(g.at(PointS(-60,-20,-20)) = 5);

	ATALWAYS(g.isInBounds(PointS(-100,0,0)) == 0,"Unexpected out-of-bounds condition");

	g.setDomainAction(Grid3D::domain_expand,Grid3D::lBound,0);
	g.setDomainAction(Grid3D::domain_subst,Grid3D::uBound,0);

	g.at(PointS(-60,-20,-20)) = 5;

	ATALWAYS(g.at(PointS(-60,-20,-20)) == 5,"Grid failed to auto-expand");

	Grid3D g1(g);

	Grid3D::SpatialDomain g_spatialDomain = g.getSpatialDomain();

	g1.resizeAndPreserve(g_spatialDomain[0]/2,g_spatialDomain[1]/2);

	for(int i=0; i < 2; ++i) {
		ATALWAYS(!blitz::all(g1.getSpatialDomain()(i) == g.getSpatialDomain()(i)),"Spatial domain did not change after resize");
	}

	Grid3D g2 = Grid::getRegionCopyS(g,g_spatialDomain[0]/2,g_spatialDomain[1]/2);

	for(int i=0; i < 2; ++i) {
		ATALWAYS(blitz::all(g1.getSpatialDomain()(i) == g2.getSpatialDomain()(i)),"Expected identical spatial domains after RegionCopy");
	}

	// build grid with spatial step 20 and logical size (50,51,51) and initial values set to 1

	Grid3D g3(PointS(20,20,20),PointL(50,51,51),true,1);

	ATALWAYS_BEGIN(blitz::all(g3.getGridArray().shape() == PointL(50, 51, 51)));
	ATDBGERR << "Unexpected internal array shape:" << g3.getGridArray().shape() << "\n";
	ATALWAYS_END;

	for(int i=0; i < 2; ++i) {
		ATALWAYS(blitz::all(g3.getSpatialDomain()(i) == Grid3D::SpatialDomain(PointS(-500., -510., -510.),PointS(500.,  510.,  510.))(i)),\
			"Unexpected spatial domain");
	}

	for(int i=0; i < 2; ++i) {
		ATALWAYS(blitz::all(g3.getLogicalDomain()(i) == Grid3D::LogicalDomain(PointL(0, 0, 0), PointL(49,  50,  50))(i)),\
			"Unexpected logical domain");
	}

}

