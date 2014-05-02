//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include <blitz/array.h>
#include <random/uniform.h>

#include "PRODDL/Geom/gdiam_simple.hpp"

#include <iostream>

#include "gtest/gtest.h"

typedef blitz::TinyVector<double,3> bz_point3;

typedef blitz::Array<bz_point3,1> bz_arr_point3;

TEST(GdiamSimpleTest, All) {

  using namespace PRODDL::Geom::Gdiam;
  using namespace std;
  
  const int N = 1000000;
  
  bz_arr_point3 points(N);
  
  ranlib::Uniform<double> rndgen;

  for( int i = 0; i < N; i++ )
  points(i) = bz_point3(rndgen.random(),rndgen.random(),rndgen.random());

  std::cout << "points[0...10]= " << endl;
  for(int i = 0; i < 10 && i < points.size(); i++)
    std::cout << points(i) << endl;

  bz_point3 left, right;

  std::cout << "Testing simp_approx_diam_pair() - locating diameter of a set of points." << std::endl;
  
  simp_approx_diam_pair(reinterpret_cast<const double*>(points.data()), N, 
			left.data(), right.data(), 
			0.0);
  
  std::cout << "left= " << left << endl << "right= " << right << std::endl;

  std::cout << "Testing simp_approx_mvbb_grid_sample() - locating min volume bounding box" << std::endl; 

  // clear for next function call
  left = 0;
  right = 0;

  bz_arr_point3 directions(3);
  bz_point3 zero_point;
  zero_point = 0.;
  directions = zero_point;

  int grid_n = 5;
  int sample_n = 400;

  simp_approx_mvbb_grid_sample(reinterpret_cast<const double*>(points.data()), N, 
			       left.data(), right.data(),
			       reinterpret_cast<double*>(directions.data()),
			       grid_n, sample_n);

  std::cout << "left= " << left << std::endl << "right= " << right << std::endl;
  std::cout << "directions= " << directions << std::endl;

}
