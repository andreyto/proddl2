//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_CTESTC_H__
# define AT_CTESTC_H__

// macro definition that runs C++ function with a predefined name 
// 'main_test().'

#include "PRODDL/Testing/exception.hpp"

#include "PRODDL/Testing/config.hpp"

#include <algorithm>
#include <iostream>
#include <string>

namespace PRODDL { namespace Testing {

  std::string testDataDir;
  std::string testOutputDir=".";

} } // namespace PRODDL { namespace Testing { 


inline char* get_simple_cmd_option(char ** begin, char ** end, const std::string & option) {
    char ** it = std::find(begin, end, option);
    if (it != end && ++it != end) {
        return *it;
    }
    return 0;
}

inline int init_testing(int argc, char **argv) {
  const char *opt = get_simple_cmd_option(argv,argv+argc,"--test-data-dir");
  if (opt) {
      PRODDL::Testing::testDataDir=opt;
  }
  //This should will be returned as a process exit code if not 0
  if(PRODDL::Testing::testDataDir.empty()) {
      std::cerr << "Usage: --test-data-dir <directory with test input data files>\n";
      return 1;
  }
  return 0;
}


#define CTEST_MODULE() \
int main_test(); \
int main(int argc, char **argv) \
{ \
    return init_testing(argc,argv) || main_test(); \
}

#endif // AT_CTESTC_H__
