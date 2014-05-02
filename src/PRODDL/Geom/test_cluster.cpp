//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
// Test geom/pairdist.hpp by comparing results of strainghtforward double-loop algorithm
// with the actual implementation.

#include "PRODDL/Geom/cluster.hpp"

#include <vector>

#include <algorithm>

#include "PRODDL/Common/math.hpp"

#include "PRODDL/Common/debug.hpp"

#include <boost/timer.hpp>

#include "PRODDL/Common/bz_vect_ext.hpp"

#include "PRODDL/Geom/uniform_points_sphere2.hpp"

#include "PRODDL/Testing/exception.hpp"

#include "gtest/gtest.h"

namespace PRODDL {


	template<typename T_num, int M_hash>
	struct TestClusterMatrix {

		typedef PRODDL::ClusterMatrix<T_num,M_hash> ClustMatr;

		typedef typename ClustMatr::Matrix Matrix;
		typedef typename ClustMatr::Floats Floats;
		typedef typename ClustMatr::Ints Ints;



		void generateRandomSphereCluster(const Floats& center, T_num radius, Matrix& m) {
			using namespace blitz;
			T_num radius2 = radius*radius;
			ranlib::UniformOpenClosed<T_num> rndgen;
			int N_rows = m.rows();
			int M_cols = m.cols();
			for(int i = 0; i < N_rows; i++) {
				while(true) {
					Floats row = m(i,Range::all());
					// generate randomly oriented vector within cube [-radius,radius]
					for(int j=0; j < M_cols; j++) row(j) = (rndgen.random() - 0.5)*2*radius;
					// if point is within a sphere of required radius, break out of the loop
					// and process the next point
					T_num r2 = blitz_ext::dotSelf(row);
					if(r2 <= radius2) {
						row += center;
						break;
					}
				}
			}
		}

		void generateRandomSphereClusters(int N_clust, T_num radius, 
			T_num minDist, T_num maxDist,
			Matrix& m) {

				using namespace blitz;
				int N_rows = m.rows();
				int M_cols = m.cols();

				// we should be able to divide matrix into even size clusters
				ATLOG_ASSERT_1( N_rows % N_clust == 0 );

				int nClusterSize = N_rows / N_clust;

				// box to accomodate N_clust with max distance
				// in M_cols dimensions
				// between centers maxDist
				T_num boxSide = maxDist * pow(T_num(N_clust),T_num(1.)/M_cols);

				T_num minDist2 = minDist * minDist;

				Matrix centers(N_clust,M_cols);

				ranlib::UniformOpenClosed<T_num> rndgen;

				for(int i_clust = 0; i_clust < N_clust; ) {

					Floats center_i = centers(i_clust,Range::all());
					// generate random point within cube [0,boxSide] x M_cols
					for(int j=0; j < M_cols; j++) 
						center_i(j) = rndgen.random()*boxSide;
					bool done_i = true;
					for(int j_clust = 0; j_clust < i_clust; j_clust++) {
						T_num r2 = blitz_ext::dotSelf(center_i - centers(j_clust,Range::all()));
						if( r2 < minDist2 ) {
							done_i = false;
							break;
						}
					}
					if( done_i ) {
						Matrix m_clust = m(Range(i_clust*nClusterSize,(i_clust+1)*nClusterSize-1),
							Range::all());
						generateRandomSphereCluster(center_i, radius, m_clust);	  
						i_clust++;
					}
				}
		}

		void testClustering(int N_rows, int M_cols) {

			using namespace blitz;

			ATDBGOUT << "Testing with " << N_rows << " x " << M_cols << " matrix, " << 
				"hashing " << M_hash << " columns.\n";

			Matrix m(N_rows,M_cols);

			Floats weights(N_rows);

			weights = 1.;

			T_num radius = 1;

			T_num minDist = radius * 5;
			T_num maxDist = radius * 7;

			T_num cutoff  = radius * 2;

			int N_clust = 10;

			ClustMatr clust;

			{
				boost::timer t;

				generateRandomSphereClusters(N_clust, radius, minDist, maxDist, m);

				ATDBGOUT << "Generated " << N_clust << " test clusters ";
				ATDBGOUT << "In " << t.elapsed() << " sec.\n";
			}

			{
				boost::timer t;

				clust.cluster(m,weights,cutoff);

				int N_clust_found = clust.numClusters();

				ATLOG_ASSERT_1(N_clust_found == N_clust);

				ATDBGOUT << "Done clustering; Found " << N_clust_found << " clusters. ";
				ATDBGOUT << "In " << t.elapsed() << " sec.\n";
			}
		}

	}; // struct TestClusterMatrix

} // namespace PRODDL


TEST(ClusterTest, All) {
	typedef float T_num;
	const int M_hash = 9;
	PRODDL::TestClusterMatrix<T_num,M_hash> tester;
	const int N_rows = 1000, M_cols = 9;
	tester.testClustering(N_rows,M_cols);
}


