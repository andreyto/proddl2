//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Common/nd_index_iter.hpp"

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/to_string.hpp"

#include "gtest/gtest.h"

TEST(NdIndexIterTest, All) {

  using namespace PRODDL;

  int arr[4][4][5];
  int lbound[3] = {1,1,1};
  int ubound[3] = {3,3,4};
  int ind[3];

  typedef index_mover<3> ind_mover;

  int i = 1;
  for(ind_mover::first(ind,lbound); ind_mover::is_good(ind,ubound); ind_mover::next(ind,lbound,ubound),i++) {
     arr[ind[0]][ind[1]][ind[2]] = i;
  }

  i = 1;
  for(ind_mover::before_first(ind,lbound); ind_mover::next(ind,lbound,ubound); i++) {
    int val = arr[ind[0]][ind[1]][ind[2]];
    ATALWAYS(val == i,"Expected value: "+to_string(i)+". Found value: "+to_string(val));
  }  
  
  i = 1;
  for(ind[0] = lbound[0]; ind[0] < ubound[0]; ind[0]++) {
    for(ind[1] = lbound[1]; ind[1] < ubound[1]; ind[1]++) {
      for(ind[2] = lbound[2]; ind[2] < ubound[2]; ind[2]++, i++) {
	int val = arr[ind[0]][ind[1]][ind[2]];
	ATALWAYS(val == i,"Expected value: "+to_string(i)+". Found value: "+to_string(val));
	ATDBGOUT << val << '\t';
      }
      ATOUTENDL();
    }
    ATOUTENDL();
  }
  std::cout << "Ok" << std::endl;

}

