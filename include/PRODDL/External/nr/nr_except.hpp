#ifndef AT_NR_EXCEPT_H__
# define AT_NR_EXCEPT_H__

/* 
   Andrei Tovchigrechko (2001)
   Exception to throw if error occured inside adapted to C++ NR routine
*/

#include <exception>
#include <string>

namespace PRODDL { namespace nr {


 class nr_error : public std::exception {
 protected:
   std::string m_msg;
 public:
   nr_error(const std::string& msg) : 
     m_msg(msg) {
   }
   virtual const char* what() const throw() {
     return m_msg.c_str();
   }
   virtual ~nr_error() throw() { }
 };

} } // namespace PRODDL::nr

#endif // AT_NR_EXCEPT_H__

