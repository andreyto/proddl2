//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Common/logger.hpp"

#include <iostream>

#include "gtest/gtest.h"

using namespace PRODDL;

namespace {


    class my_null_streambuf : public std::streambuf
    {
        public:

            my_null_streambuf()  {}
            ~my_null_streambuf() {}

      //int pubsync() { return 0; }

        protected:

            int overflow(int) { return 0; }
            int sync()        { return 0; }
    };


std::ostream my_null_ostream(new my_null_streambuf());


}

TEST(LoggerTest, All) {
}
