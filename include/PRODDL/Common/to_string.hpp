//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef COMMON_TO_STRING_H_
#define COMMON_TO_STRING_H_

/*

   Code based on boost::lexical_cast<std::string>, which is broken on gcc 2.96

*/

#include <boost/config.hpp>
# ifndef BOOST_NO_STRINGSTREAM
#  include <sstream>
# else
#  include <strstream>
# endif
#include <typeinfo>

namespace PRODDL
{
    // exception used to indicate runtime to_string_cast failure
    class bad_to_string_cast : public std::bad_cast
    {
    public:
        // constructors, destructors, and assignment operator defaulted

        // function inlined for brevity and consistency with rest of library
        virtual const char * what() const throw()
        {
            return "bad to_string cast: "
                   "source type value could not be interpreted as target";
        }
    };

    template<typename Source>
    std::string to_string(const Source& arg)
    {
# ifndef BOOST_NO_STRINGSTREAM
        std::ostringstream interpreter;
# else
        std::ostrstream interpreter; // for out-of-the-box g++ 2.96
# endif

        if(!(interpreter << arg) )
            throw bad_to_string_cast();

        return interpreter.str();
    }
} // namespace PRODDL

#endif // COMMON_TO_STRING_H_
