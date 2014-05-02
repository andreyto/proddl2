//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PIPE_STREAM_EXCEPT_H__
# define AT_PIPE_STREAM_EXCEPT_H__

// Exception class for functions from pipe_stream.hpp

#include <exception>
#include <string>

namespace PRODDL { 

  class pipe_stream_error : public std::exception {
  protected:
    std::string m_msg;
  public:
    pipe_stream_error(const std::string& msg)  throw (): 
      m_msg(msg) {
    }
    virtual const char* what() const throw (){
      return m_msg.c_str();
    }
    virtual ~pipe_stream_error() throw (){}
  };

} // namespace PRODDL

#endif // AT_PIPE_STREAM_EXCEPT__
