#include "gtest/gtest.h"

#include "PRODDL/Blitz/bztinyvec.hpp"
#include <blitz/array.h>
#include <blitz/tinymat2.h>
#include <blitz/tinymat2.cc>

#include <iostream>

using namespace blitz;
using namespace std;

TEST(BzTinyTest, ProdVec)
{
	TinyMatrix<int,3,2> A;
    A = 1;
	TinyVector<int,2> b;
	b = 1;
	TinyVector<int,3> res;
	res = 2;
	EXPECT_TRUE(all(blitz_ext::product(A,b) == res));
}

TEST(BzTinyTest, ProdMatr)
{
	TinyMatrix<int,2,3> A;
    A = 1;
	TinyMatrix<int,3,2> B;
	B = 1;
	TinyMatrix<int,2,2> res;
	res = 3;
	EXPECT_TRUE(all(blitz_ext::product(A,B) == res));
}
