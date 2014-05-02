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

// Internally, dbg:: library is controlled by the presence DBG_ENABLED macro.
// We will always define it, unless out macro DBG_DISABLED is passed from the outside,
// in which case we also set ATLOG_LEVEL to 0

#ifdef DBG_DISABLED
# undef DBG_ENABLED
#endif

#ifndef DBG_DISABLED
#ifndef DBG_ENABLED
#define DBG_ENABLED
#endif
#endif

#include "PRODDL/External/Dbg/dbg.h"
//#include "PRODDL/External/Dbg/debugbreak.h"

#include <iostream>

// Constants: threshold levels of logging output

#define ATLOG_LEVEL_0 0 // User info that is almost always on

#define ATLOG_LEVEL_1 1 // User info that is almost always on

#define ATLOG_LEVEL_2 2 // Detailed info for the user

#define ATLOG_LEVEL_3 3 // Info for the developer: tracing etc

#define ATLOG_LEVEL_4 4 // Abandand info for the developer - output from inner loops etc

#define ATLOG_LEVEL_5 5 // Abandand info for the developer - output from inner loops etc

#define ATLOG_LEVEL_6 6 // Abandand info for the developer - output from inner loops etc

// Compile time level of debugging output. ATLOG_LEVEL has precedence
// over run time level Logger::RunTimeLevel in ATLOG_SWITCH_X macros. 
// The run time level is checked only if the compile time level is satisfied.
// Thus, the programm can be compiled with higher level debugging output
// eliminated altogether, and the remaining levels can be chosen
// during run time.

#ifndef ATLOG_LEVEL
# define ATLOG_LEVEL ATLOG_LEVEL_1
#endif // ATLOG_LEVEL

#ifndef DBG_ENABLED
# undef ATLOG_LEVEL
# define ATLOG_LEVEL ATLOG_LEVEL_0
#endif


namespace PRODDL {

  class Logger {

  public:

    static void init() {

      //enable all messages for all sources
      dbg::enable_all(dbg::all, true);

      dbg::attach_ostream(dbg::all, std::cout);
      // now all 'info' messages go to cout

      //dbg::enable_level_prefix(true);

      //ATOUTLINE("Logger() done.");

      RunTimeLevel = ATLOG_LEVEL_1;

      disableExceptions();

    }

    static void disableExceptions() {

      dbg::set_assertion_behaviour(dbg::all,dbg::assertions_abort);
    
    }

    static void enableExceptions() {

      dbg::set_assertion_behaviour(dbg::all,dbg::assertions_throw);

    }

    static int getRunTimeLevel() {

      return RunTimeLevel;

    }

    static void setRunTimeLevel(int _level) {

      if(_level >= ATLOG_LEVEL_3) {

	dbg::enable_all(dbg::tracing, true);	

      }
      else {

	dbg::enable_all(dbg::tracing, false);

      }

      RunTimeLevel = _level;

    }

  protected:

    // Sets the run time level of logging output
    static int RunTimeLevel;

  };

} // namespace PRODDL


#ifdef DBG_ENABLED

#   define ATLOGVAR(var) ""#var" = " << (var) << '\t'

#else

#   define ATLOGVAR(var) (0)

#endif // DBG_ENABLED





#if ( ATLOG_LEVEL >= ATLOG_LEVEL_1 ) 

#   define ATLOG_SWITCH_1(x)  if( PRODDL::Logger::getRunTimeLevel() >= ATLOG_LEVEL_1 ) {x;}

#   define ATTRACE_SWITCH_1(x) x

#else

#   define ATLOG_SWITCH_1(x)  (void(0))

#   define ATTRACE_SWITCH_1(x) (void(0))

#endif // ( ATLOG_LEVEL >= ATLOG_LEVEL_1 )



#if ( ATLOG_LEVEL >= ATLOG_LEVEL_2 ) 

#   define ATLOG_SWITCH_2(x)   if( PRODDL::Logger::getRunTimeLevel() >= ATLOG_LEVEL_2 ) {x;}

#   define ATTRACE_SWITCH_2(x) x

#else

#   define ATLOG_SWITCH_2(x)  (void(0))

#   define ATTRACE_SWITCH_2(x) (void(0))

#endif // ( ATLOG_LEVEL >= ATLOG_LEVEL_2 )



#if ( ATLOG_LEVEL >= ATLOG_LEVEL_3 ) 

#   define ATLOG_SWITCH_3(x)   if( PRODDL::Logger::getRunTimeLevel() >= ATLOG_LEVEL_3 ) {x;}

#   define ATTRACE_SWITCH_3(x) x

#else

#   define ATLOG_SWITCH_3(x)  (void(0))

#   define ATTRACE_SWITCH_3(x) (void(0))

#endif // ( ATLOG_LEVEL >= ATLOG_LEVEL_3 )


#if ( ATLOG_LEVEL >= ATLOG_LEVEL_4 ) 

#   define ATLOG_SWITCH_4(x)   if( PRODDL::Logger::getRunTimeLevel() >= ATLOG_LEVEL_4 ) {x;}

#   define ATTRACE_SWITCH_4(x) x

#else

#   define ATLOG_SWITCH_4(x)  (void(0))

#   define ATTRACE_SWITCH_4(x) (void(0))

#endif // ( ATLOG_LEVEL >= ATLOG_LEVEL_4 )


#if ( ATLOG_LEVEL >= ATLOG_LEVEL_5 ) 

#   define ATLOG_SWITCH_5(x)   if( PRODDL::Logger::getRunTimeLevel() >= ATLOG_LEVEL_5 ) {x;}

#   define ATTRACE_SWITCH_5(x) x

#else

#   define ATLOG_SWITCH_5(x)  (void(0))

#   define ATTRACE_SWITCH_5(x) (void(0))

#endif // ( ATLOG_LEVEL >= ATLOG_LEVEL_5 )


#if ( ATLOG_LEVEL >= ATLOG_LEVEL_6 ) 

#   define ATLOG_SWITCH_6(x)   if( PRODDL::Logger::getRunTimeLevel() >= ATLOG_LEVEL_6 ) {x;}

#   define ATTRACE_SWITCH_6(x) x

#else

#   define ATLOG_SWITCH_6(x)  (void(0))

#   define ATTRACE_SWITCH_6(x) (void(0))

#endif // ( ATLOG_LEVEL >= ATLOG_LEVEL_6 )



#define ATLOG_ASSERT_1(x) ATLOG_SWITCH_1(dbg::assertion(dbg::error,DBG_ASSERTION(x)))
#define ATLOG_ASSERT_2(x) ATLOG_SWITCH_2(dbg::assertion(dbg::error,DBG_ASSERTION(x)))
#define ATLOG_ASSERT_3(x) ATLOG_SWITCH_3(dbg::assertion(dbg::error,DBG_ASSERTION(x)))
#define ATLOG_ASSERT_4(x) ATLOG_SWITCH_4(dbg::assertion(dbg::error,DBG_ASSERTION(x)))
#define ATLOG_ASSERT_5(x) ATLOG_SWITCH_5(dbg::assertion(dbg::error,DBG_ASSERTION(x)))
#define ATLOG_ASSERT_6(x) ATLOG_SWITCH_6(dbg::assertion(dbg::error,DBG_ASSERTION(x)))



#define ATLOG_OUT_1(x) ATLOG_SWITCH_1(dbg::out(dbg::info) << dbg::indent() << x << " At: " << DBG_HERE << '\n')
#define ATLOG_OUT_2(x) ATLOG_SWITCH_2(dbg::out(dbg::info) << dbg::indent() << x << " At: " << DBG_HERE << '\n')
#define ATLOG_OUT_3(x) ATLOG_SWITCH_3(dbg::out(dbg::info) << dbg::indent() << x << " At: " << DBG_HERE << '\n')
#define ATLOG_OUT_4(x) ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() << x << " At: " << DBG_HERE << '\n')
#define ATLOG_OUT_5(x) ATLOG_SWITCH_5(dbg::out(dbg::info) << dbg::indent() << x << " At: " << DBG_HERE << '\n')
#define ATLOG_OUT_6(x) ATLOG_SWITCH_6(dbg::out(dbg::info) << dbg::indent() << x << " At: " << DBG_HERE << '\n')

// ATLOG_TRACE_OBJECT is not a macro, it is a literal variable name that will be created
#define ATLOG_TRACE_1 ATTRACE_SWITCH_1(dbg::trace ATLOG_TRACE_OBJECT(DBG_HERE))
#define ATLOG_TRACE_2 ATTRACE_SWITCH_2(dbg::trace ATLOG_TRACE_OBJECT(DBG_HERE))
#define ATLOG_TRACE_3 ATTRACE_SWITCH_3(dbg::trace ATLOG_TRACE_OBJECT(DBG_HERE))
#define ATLOG_TRACE_4 ATTRACE_SWITCH_4(dbg::trace ATLOG_TRACE_OBJECT(DBG_HERE))
#define ATLOG_TRACE_5 ATTRACE_SWITCH_5(dbg::trace ATLOG_TRACE_OBJECT(DBG_HERE))
#define ATLOG_TRACE_6 ATTRACE_SWITCH_6(dbg::trace ATLOG_TRACE_OBJECT(DBG_HERE))

#if ( ATLOG_LEVEL >= ATLOG_LEVEL_1 )

#   define ATLOG_STD_EXCEPTIONS_TRY()  try { ATLOG_OUT_4("Entered try statement for standard exceptions");

#   define ATLOG_STD_EXCEPTIONS_CATCH()  \
	ATLOG_OUT_4("End of try statement for standard exceptions"); \
	} \
	catch (const std::exception& x) { \
	  dbg::out(dbg::error) << dbg::indent() << "Std::exception: " << ATLOGVAR(x.what()) << "\n"; \
		throw; \
	}	\
    catch (...) { \
	  dbg::out(dbg::error) << dbg::indent() << "Unknown exception caught.\n"; \
		throw; \
    }

#else

#   define ATLOG_STD_EXCEPTIONS_TRY()  {}

#   define ATLOG_STD_EXCEPTIONS_CATCH()  {}

#endif // ( ATLOG_LEVEL >= ATLOG_LEVEL_1 )

#endif // PRODDL_LOGGER_H__
