//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_COMMON_DEBUG_H__
#define AT_COMMON_DEBUG_H__

/*
   Copyright: Andrei Tovchigrechko, 2000

   Some basic debugging macros extending standard 'assert()' functionality.

*/

#include <iostream> //for debug ouput
#include <exception>

#include <cstdlib>

#if defined(ATDEBUG)
#include "PRODDL/External/Dbg/debugbreak.h"
#endif

// for printing ranges to standard output
#include <iterator>
#include <algorithm>

#ifndef ATDEBUG_LEVEL
#define ATDEBUG_LEVEL 6
#endif

#define ATDBGOUT (std::cout)
#define ATDBGERR (std::cerr)

struct dbg_exception : public std::exception {

};

#if defined(ATDEBUG) && defined(ATDEBUG_BREAK)
# define AT_THROW(_ex) (debug_break(),throw (_ex))
#else
# define AT_THROW(_ex) throw (_ex)
#endif

// inline void ATDBREAK(std::ostream& os) {
//   //((char*)0)[0] = 1;
//     //std::exit(1);
//     throw dbg_exception();
// }

#define ATDBREAK(_out_stream) (AT_THROW(dbg_exception()))


#if defined(ATDEBUG)
#define ATASSERT(st,msg) \
	(!(st)? \
		(ATDBGERR<<"assertion failed: LINE "<<__LINE__<<"; FILE "<<__FILE__\
	<<"; statement: "<<""#st""<<"; Message: "<<(msg)<<'\n',ATDBREAK(ATDBGERR)):(void)0)

#else
#define ATASSERT(st,msg) ((void)0)   
#endif

#define ATALWAYS(st,msg) \
	(!(st)? \
		(ATDBGERR<<"assertion failed: LINE "<<__LINE__<<"; FILE "<<__FILE__\
	<<"; statement: "<<""#st""<<"; Message: "<<(msg)<<'\n',ATDBREAK(ATDBGERR)):(void)0)

// if statement 'st' is false, print the message 'msg' and code location to ATDBGERR and
// throw exception '_exception_type', whose ctor must take string 'msg' as an argument
#define ATALWAYS_EX(st,msg,_exception_type) \
	(!(st)? \
		(ATDBGERR<<"assertion failed: LINE "<<__LINE__<<"; FILE "<<__FILE__\
	<<"; statement: "<<""#st""<<"; Message: "<<(msg)<<'\n',AT_THROW(_exception_type(msg))):(void)0)


#define ATALWAYS_BEGIN(st__) if( ! (st__) ) { \
        ATDBGERR <<"assertion failed: LINE "<<__LINE__<<"; FILE "<<__FILE__\
	<<"; statement: "<<""#st__""<<"; Messages:" << '\n'; 


#define ATALWAYS_END ATDBGERR << "End of assertion block" << '\n'; ATDBREAK(ATDBGERR); }

#if defined(ATDEBUG)
# define ATASSERT_BEGIN(st__) if( ! (st__) ) { \
        ATDBGERR <<"assertion failed: LINE "<<__LINE__<<"; FILE "<<__FILE__\
	<<"; statement: "<<""#st__""<<"; Messages:" << '\n'; 

# define ATASSERT_END ATDBGERR << "End of assertion block" << '\n'; ATDBREAK(ATDBGERR); }
#else
# define ATASSERT_BEGIN(st__) if(0) {
# define ATASSERT_END }
#endif

#define ATDIE(msg) \
	(ATDBGERR<<"requirement failed: LINE "<<__LINE__<<"; FILE "<<__FILE__\
	<<"; Message: "<<(msg)<<'\n',ATDBREAK(ATDBGERR))


#if defined(ATDEBUG)

#define ATDBGVAR(var) (ATDBGOUT<<""#var" = "<<(var)<<" in line: "<<__LINE__\
	<<" in file: "<<__FILE__<<'\n')

#else

#define ATDBGVAR(var) ((void)0)

#endif


#define ATOUTVAR(var) (ATDBGOUT<<""#var" = "<<(var)<<'\t')
#define ATOUTRANGE(_at__start,_at__end,_at__type) ((ATDBGOUT<<"["#_at__start","#_at__end") = "),\
         (std::copy((_at__start),(_at__end),std::ostream_iterator<_at__type>(ATDBGOUT, "\t"))))
#define ATOUTENDL()     (ATDBGOUT<<'\n')

#define ATSOUTVAR(var,out_stream__) ((out_stream__)<<""#var" = "<<(var)<<'\t')
#define ATSOUTENDL(out_stream__)     ((out_stream__)<<'\n')

#if defined(ATDEBUG)

#define ATOUTLINE(var) (ATDBGOUT<<(var)<<" in line: "<<__LINE__\
	<<" in file: "<<__FILE__<<'\n')

#else

#define ATOUTLINE(var) ((void)0)

#endif


#if defined(ATDEBUG)

#define ATPROGRESS(i) (ATDBGOUT<<"PROGRESS "<<(i)<<" in line: "<<__LINE__\
	<<" in file: "<<__FILE__<<'\n')

#else

#define ATPROGRESS(i) ((void)0)

#endif


//You can wrap any statement by TRACE
//macros. If ATDEBUG is not defined, statement
//will be just evaluated. Otherwise, it
//will be evaluated and it's definition (not
//it's result) will be printed.

#if defined(ATDEBUG)

#define ATTRACE(st) \
((ATDBGOUT<<"TRACE: LINE "<<__LINE__<<";FILE "<<__FILE__\
	<<"; statement "<<#st<<'\n'),(st),(void)0)

#else

#define ATTRACE(st) ((st),(void)0)

#endif


#if defined(ATDEBUG)

  template<class Cont, class Ret, class Ind>
  Ret& 
  AT_TAKE_SQUARE(Cont& c, const Ind& i, const char *file, int line,const Ret&) {
    if( i >= c.size() || i < 0 ) {
      char *p = 0;
      *p  = '\0';
      ATDBGOUT << "Index out of range at " << file << " : " << 
	line << " : " << c.size() << " : " << i << "\n";
      throw dbg_exception();
    }
    return c[i];
  }

  template<class Cont, class Ret, class Ind>
  const Ret& 
  AT_TAKE_SQUARE(const Cont& c, const Ind& i, const char *file, int line, const Ret&) {
    if( i >= c.size() || i < 0 ) {
      ATDBGOUT << "Index out of range at " << file << " : " << 
	line << " : " << c.size() << " : " << i << "\n";
      AT_THROW(dbg_exception());
    }
    return c[i];
  }

  template<class Cont, class Ret, class Ind>
  Ret& 
  AT_TAKE_ROUND(Cont& c, const Ind& i, const char *file, int line,const Ret&) {
    if( i >= c.size() || i < 0 ) {
      ATDBGOUT << "Index out of range at " << file << " : " << 
	line << " : " << c.size() << " : " << i << "\n";
      AT_THROW(dbg_exception());
    }
    return c(i);
  }

  template<class Cont, class Ret, class Ind>
  Ret
  AT_TAKE_ROUND(const Cont& c, const Ind& i, const char *file, int line, const Ret&) {
    if( i >= c.size() || i < 0 ) {
      ATDBGOUT << "Index out of range at " << file << " : " << 
	line << " : " << c.size() << " : " << i << "\n";
      AT_THROW(dbg_exception());
    }
    return c(i);
  }


#define ATTAKE_DBG(_cont,_ind) AT_TAKE_SQUARE(_cont,_ind,__FILE__,__LINE__,(_cont)[0])

#define ATTAKE_DBG_R(_cont,_ind) AT_TAKE_ROUND(_cont,_ind,__FILE__,__LINE__,(_cont)(0))

#else

#define ATTAKE_DBG(_cont,_ind) (_cont)[_ind]

#define ATTAKE_DBG_R(_cont,_ind) (_cont)(_ind)

#endif


#endif // AT_COMMON_DEBUG_H__

