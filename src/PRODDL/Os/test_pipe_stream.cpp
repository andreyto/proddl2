//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Os/pipe_stream.hpp"

#include <string>

#include "PRODDL/Testing/exception.hpp"

#include "gtest/gtest.h"

/// @todo make the test portable across OSes

TEST(PipeStreamTest, All) {
  PRODDL::opipe_stream out_pipe("sort");
  for(int i = 9; i >=0; i--)
    out_pipe << i << std::endl;
  PRODDL::ipipe_stream in_pipe("echo 'abdication' | tr 'abcd' 'wxyz'");
  std::string expected = "wxziywtion";
  std::string received;
  in_pipe >> received;
  if( received != expected ) {
    throw PRODDL::Testing::test_error("ipipe_stream input does not match. Expected: " + 
				      expected + ". Received: "+received);
  }
}
