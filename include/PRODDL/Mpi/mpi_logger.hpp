//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_MPI_LOGGER_H__
#define PRODDL_MPI_LOGGER_H__

/*

*/

#include "PRODDL/Common/logger.hpp"

#include "PRODDL/Mpi/mpi.hpp"

#include <sstream>

namespace PRODDL {

  class MpiLogger : public Logger {

  public:


    // This method makes MPI calls, and so must be called AFTER the call to MPI_Init()

    static
    void
    postInit() {

      std::ostringstream out;

      out << "R"
	  << std::setw(3) << std::setprecision(3) 
	  << std::setfill('0') << OOMPI_COMM_WORLD.Rank() << ' ' << std::ends;

      // set prefix to the Rank of this process
      dbg::set_prefix(out.str().c_str());

    }

  };

} // namespace PRODDL

#endif // PRODDL_MPI_LOGGER_H__
