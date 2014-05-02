//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Common/queue.hpp"

#include "PRODDL/Common/debug.hpp"

#include <iterator>

#include <algorithm>

#include "gtest/gtest.h"

using namespace PRODDL;

TEST(BoundPriorityQueueTest, All) {

  BoundPriorityQueue<int> q;

  q.init(4,10);

  q.push(1);

  q.push(5);

  q.push(11);

  q.push(12);

  q.push(4);

  q.push(3);

  q.push(-5);

  q.push(2);

  const std::vector<int>& data = q.data();

  for(int i = 0; i < data.size(); i++) {

    ATOUTVAR(i); ATOUTVAR(data[i]);

  }

  ATOUTENDL();

  ATALWAYS(q.size() == 4,"");

  ATALWAYS(q.top() == 3,"");

  q.sort();

  for(int i = 0; i < data.size(); i++) {

    ATOUTVAR(i); ATOUTVAR(data[i]);

  }

  ATALWAYS(data[0] == -5,"");
  ATALWAYS(data[1] == 1,"");
  ATALWAYS(data[2] == 2,"");
  ATALWAYS(data[3] == 3,"");

  try {

    q.push(2);

  }
  catch(BoundPriorityQueueError msg) {

    ATOUTVAR(msg.what()); ATOUTENDL();

  }

  q.clear();

  ATALWAYS(q.size() == 0,"");

  q.push(2);

  ATALWAYS(q.top() == 2,"");

  q.pop();

  ATALWAYS(q.size() == 0,"");

  int a[5] = { 5, 4, 12, 1, 2 };

  std::copy(a,a+5,std::back_inserter(q));

  ATALWAYS(q.size() == 4,"");

}

