//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_LOGGER_H__
#define PRODDL_LOGGER_H__

/*

*/

#include <iostream>

namespace PRODDL {

  class Logger {

  public:

    Logger() {

      logLevel = 0;

      lastLogLevel = 1;

      out.attach(std::cout);

    }

    Logger& operator()(int _logLevel) {

      lastLogLevel = _logLevel;

      return *this;

    }

    template<class T>
    Logger& operator<< (const T& x) {

      if( lastLogLevel <= logLevel ) {

	out << x;

      }

      return *this;

    }

    int level() const {

      return logLevel;

    }

    void level(int newLevel) {

      logLevel = newLevel;

    }

  protected:

    int logLevel;

    int lastLogLevel;

    ostream out;

  };

} // namespace PRODDL

#endif // PRODDL_LOGGER_H__
