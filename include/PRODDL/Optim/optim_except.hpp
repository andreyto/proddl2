//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_OPTIM_EXCEPT_H__
# define AT_OPTIM_EXCEPT_H__

// Exception classes for objects in geom namespace

#include <exception>

namespace PRODDL { namespace Optim {

  class optim_error : public std::exception {
  protected:
    std::string m_msg;
  public:
    optim_error(const std::string& msg)  throw (): 
      m_msg(msg) {
    }
    virtual ~optim_error() throw (){}
    virtual const char* what() const  throw (){
      return m_msg.c_str();
    }
  };

  class optim_param_error : public optim_error {
  public:
    optim_param_error(const std::string& msg)  throw (): 
      optim_error(msg) {
    }
    virtual ~optim_param_error() throw (){}
  };

} } // namespace PRODDL::Optim

#endif // AT_OPTIM_EXCEPT_H__
