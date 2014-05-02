//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_TESTING_EXCEPTION_H__
#define AT_TESTING_EXCEPTION_H__

#include <exception>

#include <string>

namespace PRODDL { namespace Testing {

  class test_error : public std::exception {
  protected:
    std::string m_msg;
  public:
    test_error(const std::string& msg) : 
      m_msg(msg) {
    }
    virtual const char* what() const throw() {
      return m_msg.c_str();
    }
    virtual ~test_error() throw() { }
  };

} } // namespace PRODDL::Testing


#endif // AT_TESTING_EXCEPTION_H__
                                                                                                     
