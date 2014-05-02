//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_EXCEPT_H__
# define PRODDL_EXCEPT_H__

// Exception classes for objects in PRODDL namespace

#include <exception>
#include <string>


#define PRODDL_EXCEPT_DECL(_except_name) \
  class _except_name : public docktk_error { \
  public: \
    _except_name(const std::string& msg)  throw (): \
      docktk_error(msg) { \
    } \
    virtual ~_except_name() throw (){} \
  }


namespace PRODDL {

  class docktk_error : public std::exception {
  protected:
    std::string m_msg;
  public:
    docktk_error(const std::string& msg)  throw (): 
      m_msg(msg) {
    }
    virtual const char* what() const throw (){
      return m_msg.c_str();
    }
    virtual ~docktk_error() throw (){}
  };

  class potential_error : public docktk_error {
  public:
    potential_error(const std::string& msg)  throw (): 
      docktk_error(msg) {
    }
    virtual ~potential_error() throw (){}
  };

  class tolerance_error : public docktk_error {
  public:
    explicit tolerance_error(const std::string& msg):
      docktk_error(msg)
    {}
  };

  class not_supported_error : public docktk_error {
  public:
    not_supported_error(const std::string& msg)  throw (): 
      docktk_error(msg) {
    }
    virtual ~not_supported_error() throw (){}
  };

  PRODDL_EXCEPT_DECL(io_error);

  PRODDL_EXCEPT_DECL(argument_error);

  PRODDL_EXCEPT_DECL(size_error);
  
  PRODDL_EXCEPT_DECL(workflow_error);

} // namespace PRODDL

#endif // PRODDL_EXCEPT_H__
