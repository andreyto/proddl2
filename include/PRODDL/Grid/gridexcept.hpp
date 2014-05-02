//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef ATGRID_EXCEPTIONS_H__
#define ATGRID_EXCEPTIONS_H__

#include <exception>

namespace PRODDL { namespace Grid {

 class grid_error : public std::exception {

 protected:
   std::string m_msg;

 public:

   grid_error(const std::string& msg) : 
     m_msg(msg) {
   }

   virtual const char* what() const throw() {
     return m_msg.c_str();
   }

   virtual ~grid_error() throw() { }
 };

 class grid_domain_error : public grid_error {
   
 public:
   grid_domain_error(const std::string& msg) :
     grid_error(msg)
     {}
   virtual ~grid_domain_error() throw() { }
 };

 class grid_sanity_error : public grid_error {

 public:
   grid_sanity_error(const std::string& msg) :
     grid_error(msg)
     {}
   virtual ~grid_sanity_error() throw() { }
 };

} } // namespace PRODDL::Grid


#endif // ATGRID_EXCEPTIONS_H__
