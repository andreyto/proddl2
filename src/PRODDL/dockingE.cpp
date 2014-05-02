//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/docking.hpp"

#include "PRODDL/Common/logger.hpp"


namespace PRODDL {

  Options gOptions;

  Logger gLogger;

} // namespace PRODDL



namespace {

  void
  setOptions(const PRODDL::Options& options) {

    int runTimeLogLevel;

    options.getdefault("logLevel",runTimeLogLevel,ATLOG_LEVEL);

    PRODDL::Logger::setRunTimeLevel(runTimeLogLevel);

    PRODDL::gOptions = options;

  }


  typedef double T_num;

  int numSize() {

    return sizeof(T_num);
  }
  

}


