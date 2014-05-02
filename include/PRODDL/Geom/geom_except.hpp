//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_GEOM_EXCEPT_H__
# define AT_GEOM_EXCEPT_H__

// Exception classes for objects in geom namespace

#include "PRODDL/exceptions.hpp"
#include <string>

namespace PRODDL { namespace Geom {

  class geom_error : public std::exception {
  protected:
    std::string m_msg;
  public:
    geom_error(const std::string& msg)  throw (): 
      m_msg(msg) {
    }
    virtual const char* what() const throw (){
      return m_msg.c_str();
    }
    virtual ~geom_error() throw (){}
  };

  class sphere_points_error : public geom_error {
  public:
    sphere_points_error(const std::string& msg)  throw (): 
      geom_error(msg) {
    }
    virtual ~sphere_points_error() throw (){}
  };

  template<class T> class sphere_points_value_error : public geom_error {
  protected:
    const std::string name;
    const T val;
  public:
    sphere_points_value_error(const std::string& msg,const std::string& _name,const T& _val)  throw ():
      geom_error(msg),
      name(_name),
      val(_val) {
    }
    virtual ~sphere_points_value_error() throw (){}
    const std::string& getName() const {
      return name;
    }
    const T& getVal() const {
      return val;
    }
  };

  template<class T> class sphere_points_param_error : public sphere_points_value_error<int> {
  public:
    sphere_points_param_error(const std::string& msg, const std::string& _name, const T& _val)  throw (): 
      sphere_points_value_error<int>(msg,_name,_val) {
    }
    virtual ~sphere_points_param_error() throw (){}
  };

  class sphere_points_loops_error : public sphere_points_value_error<int> {
  public:
    sphere_points_loops_error(const std::string& msg, const std::string& _name, int _n_loops)  throw (): 
      sphere_points_value_error<int>(msg,_name,_n_loops) {
    }
    virtual ~sphere_points_loops_error() throw (){}
  };

  // Exceptions from code in 'geom::align' namespace

  class quatFit_error : public geom_error {
  public:
    quatFit_error(const std::string& msg)  throw (): 
      geom_error(msg) {
    }
    virtual ~quatFit_error() throw (){}
  };

  // Exception for 'transforms'

  class transforms_error : public geom_error {
  public:
    transforms_error(const std::string& msg)  throw (): 
      geom_error(msg) {
    }
    virtual ~transforms_error() throw (){}
  };

  class transforms_size_error : public transforms_error {
  public:
    transforms_size_error(const std::string& msg)  throw (): 
      transforms_error(msg) {
    }
    virtual ~transforms_size_error() throw (){}
  };

} } // namespace PRODDL::Geom

#endif // AT_GEOM_EXCEPT_H__
