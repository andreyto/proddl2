//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/Testing/exception.hpp"

#include "PRODDL/Common/debug.hpp"

#include <blitz/array.h>

#include "PRODDL/Blitz/bzutil.hpp"

#include "gtest/gtest.h"

TEST(BzExtTest, All) {

    typedef double T_num_large;
    typedef float T_num_small;
    typedef blitz::TinyVector<T_num_large,3> PointLarge;
    typedef blitz::TinyVector<T_num_small,3> PointSmall;

    PointLarge pLarge(1,2,3);
    PointSmall pSmall(4,5,6);
    pSmall = pLarge;
    ATOUTVAR(pLarge); ATOUTVAR(pSmall); ATOUTENDL();
    
    typedef blitz::Array<PointLarge,1> ArrayPointLarge;
    typedef blitz::Array<PointSmall,1> ArrayPointSmall;

    ArrayPointLarge aLarge(5);
    aLarge = 1;
    //ArrayPointSmall aSmall(blitz::cast<PointSmall>(aLarge)); //won't compile
    ArrayPointSmall aSmall(5);
    std::copy(aLarge.begin(),aLarge.end(),aSmall.begin());
    ATOUTVAR(aLarge); ATOUTVAR(aSmall); ATOUTENDL();
    ArrayPointSmall aSmall_1(blitz_ext::makeArrayWithSameStructure<PointSmall>(aLarge));
    std::copy(aLarge.begin(),aLarge.end(),aSmall_1.begin());
    ATOUTVAR(aLarge); ATOUTVAR(aSmall_1); ATOUTENDL();
    ArrayPointSmall aSmall_2(blitz_ext::getContiguousRefOrCopy<PointSmall>(aLarge));
    ATOUTVAR(aLarge); ATOUTVAR(aSmall_2); ATOUTENDL();
    ArrayPointLarge aLarge_1(blitz_ext::getContiguousRefOrCopy<PointLarge>(aLarge));
    ATOUTVAR(aLarge); ATOUTVAR(aLarge_1); ATOUTENDL();
  }

