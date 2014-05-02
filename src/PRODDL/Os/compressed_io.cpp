//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Os/compressed_io.hpp"

#include <fstream>

#include "PRODDL/External/Gzstream/gzstream.hpp"

#include <string>

#include "PRODDL/Os/pipe_stream.hpp"

namespace PRODDL {


  namespace CompressedIO {

    // types of compressed files. Unknown means any type except other
    // members of FileType enum.
    
    enum FileType { Unknown, Z, gz };

    FileType detectFileType(const char* fileName) {
      std::string name(fileName);
      std::string ext = ".Z";
      if( name.substr(name.size() - ext.size()) == ext )
	return Z;
      ext = ".gz";
      if( name.substr(name.size() - ext.size()) == ext )
	return gz;
      return Unknown;
    }

  } // namespace CompressedIO

  boost::shared_ptr<std::ostream> openCompressedOutput(const char* fileName) {

    CompressedIO::FileType file_type = CompressedIO::detectFileType(fileName);
    
    std::ostream * ret_val = 0;

    if( file_type == CompressedIO::Z ) {
      // create a pipe to "compress" program and write into it
      std::string fname(fileName);
      ret_val = new opipe_stream(("compress -c > " + fname).c_str());
    }
    else if( file_type == CompressedIO::gz ) {
      ret_val = new gzstream::ogzstream(fileName);
    }
    else {
      ret_val = new std::ofstream(fileName);
    }

    return boost::shared_ptr<std::ostream>(ret_val);

  }

  boost::shared_ptr<std::istream> openUncompressedInput(const char* fileName) {

    CompressedIO::FileType file_type = CompressedIO::detectFileType(fileName);
    
    std::istream * ret_val = 0;

    if( file_type == CompressedIO::Z ) {
     
      // create a pipe to "compress" program and read from it
      std::string fname(fileName);
      ret_val = new ipipe_stream(("uncompress -c " + fname).c_str());
    }
    else if( file_type == CompressedIO::gz ) {
      ret_val = new gzstream::igzstream(fileName);
    }
    else {
      ret_val = new std::ifstream(fileName);
    }

    return boost::shared_ptr<std::istream>(ret_val);

  }

} // namespace PRODDL
