//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_COMPRESSED_IO__
#define PRODDL_COMPRESSED_IO__

// Functions that return stream objects with on-the-fly uncompression/compression of
// most popular compression formats: gzip (.gz, using zlib) and compress (.Z, using pipes
// to/from gzip/gunzip executable).
// File type is determined by name extension (e.g. *.gz). For unknown file types the
// stream that deals with file data without any conversion is returned.

#include <iostream>

#include <boost/shared_ptr.hpp>

//#include "PRODDL/Os/compressed_io_except.hpp"

namespace PRODDL {

	boost::shared_ptr<std::ostream> openCompressedOutput(const char* fileName);

  boost::shared_ptr<std::istream> openUncompressedInput(const char* fileName);

} // namespace PRODDL

#endif // PRODDL_COMPRESSED_IO__

