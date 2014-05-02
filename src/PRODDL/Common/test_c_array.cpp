//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include <iostream>

#include "PRODDL/Common/c_array.hpp"

#include "gtest/gtest.h"

using namespace std;
using namespace PRODDL;

void f(float a[3][3]) {
}


typedef c_array_traits<float,2,3>::A_type AF3;

template<int rank,int n2, typename T> void f2(c_array_traits<T,rank,n2> x) {
}

TEST(CArrayTest, All) {

  float a[3*3];
  c_array_traits<float,2,3>::A_type a3x3 = c_array_cast<3>(a);
  for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
	{
	  a3x3[i][j] = 10*i + j;
	  cout << a3x3[i][j] << '\t';
	}
      cout << endl;
    }
  cout << endl;
  f(a3x3);
  float a1[4*3*2];
  c_array_traits<float,3,3,2>::A_type a4x3x2 = c_array_cast<3,2>(a1);
  for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < 3; j++)
	{
	  for(int k = 0; k < 2; k++)
	    {
	      a4x3x2[i][j][k] = 100*i + 10*j + k;
	      cout << a4x3x2[i][j][k] << '\t';
	    }
	  cout << endl;
	}
      cout << endl;
    }
    	
}

    
