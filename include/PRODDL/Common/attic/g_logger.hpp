//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_G_LOGGER_H__
#define PRODDL_G_LOGGER_H__

#include "PRODDL/Common/logger.hpp"

namespace PRODDL {

  extern Logger Log;

} // namespace PRODDL

#define ATLOGVAR(__var) ((Log) << ""#__var" = "<<(__var)<<'\t')

#define ATLOGENDL()     ((Log) << '\n')

#endif // PRODDL_G_LOGGER_H__
