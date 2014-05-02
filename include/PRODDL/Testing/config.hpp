//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_TESTING_CONFIG_H__
#define PRODDL_TESTING_CONFIG_H__

// Configuration data for testing procedures:
// default location for data files and temporary files

#include <string>

namespace PRODDL { namespace Testing {

  extern std::string testDataDir;
  extern std::string testOutputDir;


} } // namespace PRODDL { namespace Testing { 

#endif // PRODDL_TESTING_CONFIG_H__
