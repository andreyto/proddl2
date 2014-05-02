#include "PRODDL/Testing/ctestc.hpp"
#include "gtest/gtest.h"

GTEST_API_ int main(int argc, char **argv) {
  //this will remove all options that it parses
  ::testing::InitGoogleTest(&argc, argv);
  return init_testing(argc,argv) || RUN_ALL_TESTS();
}
