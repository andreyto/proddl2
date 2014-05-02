//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_COMPRESSED_IO_EXCEPT_H__
# define AT_COMPRESSED_IO_EXCEPT_H__

// Exception class for functions from compressed_io.hpp

#include <exception>
#include <string>

namespace PRODDL { 

  class compressed_io_error : public std::exception {
  protected:
    std::string m_msg;
  public:
    compressed_io_error(const std::string& msg)  throw (): 
      m_msg(msg) {
    }
    virtual const char* what() const throw (){
      return m_msg.c_str();
    }
    virtual ~compressed_io_error() throw (){}
  };

} // namespace PRODDL

#endif // AT_COMPRESSED_IO_EXCEPT__
