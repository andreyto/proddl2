//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_MATH_EXCEPT_H__
# define AT_MATH_EXCEPT_H__

// Exception classes for objects in Math namespace

#include <exception>
#include <string>

namespace PRODDL { namespace Math {

  class math_error : public std::exception {
  protected:
    std::string m_msg;
  public:
    math_error(const std::string& msg)  throw (): 
      m_msg(msg) {
    }
    virtual const char* what() const throw (){
      return m_msg.c_str();
    }
    virtual ~math_error() throw (){}
  };

  class fft_size_error : public math_error {
  public:
    fft_size_error(const std::string& msg)  throw (): 
      math_error(msg) {
    }
    virtual ~fft_size_error() throw (){}
  };

} } // namespace PRODDL::Math

#endif // AT_MATH_EXCEPT_H__
