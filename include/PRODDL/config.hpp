//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_CONFIG_H__
#define PRODDL_CONFIG_H__

#include <map>

#include <string>

#include <sstream>

#include "PRODDL/Common/function.hpp"


namespace PRODDL {

  // Class that provides global parameters of a running program.
  // I has a template member method that extracts a parameter with 
  // a given name into a variable of a given type.

  class Parameters {

  protected:

    typedef std::map<std::string,std::string> Dict;

  public:
    
    template<typename T_val>
    void getValue(const std::string& name, T_val& val) const {

      std::ostringstream in(getMapValue(dict,name));
      in >> val;

    }


    void getValue(const std::string& name, std::string& val) const {

      val = getMapValue(dict,name);

    }

    template<typename T_val>
    void getValue(const char* name, T_val& val) const {

      getValue(std::string(name),val);

    }


    Parameters() {

      //TODO: implement loading the string representations of parameters
      //from JSON

      dict["Params.Root"] = "/home/andrey/science/python/data";
      dict["Params.AngleGridDir"] = "/usr/data/common/var/gramm";

    }


  protected:
    
    Dict dict;

  }; // class Parameters



} // namespace PRODDL

#endif // PRODDL_CONFIG_H__
