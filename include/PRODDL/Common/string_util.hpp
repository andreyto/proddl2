//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_STRING_UTIL_H__
#define AT_STRING_UTIL_H__

// Some utility functions dealing with strings: trancation, check for 'all spaces' strings etc.

#include <string>
#include <algorithm>
#include <functional>

#include <cstring>
#include <cctype>

namespace PRODDL {

	/// return a copy of string 's' with char 'from' replaced with char 'to'
	inline std::string replace_char(const std::string& s, char from, char to) {
		std::string z(s);
		//comes from <algorithm>
		std::replace(z.begin(),z.end(),from,to);
		return z;
	}



	//define this function to use with 'std::ptr_fun' since standard
	//isspace() from string.h is often a macros

	inline bool is_space(char c)
	{
		return std::isspace(c);
	}

	// Return string with left and right spaces trancated (spaces are defined by 'isspace' from <cctype>)

	std::string strtranc(const std::string& s);

	//  char *strtranc(char *dest,const char *src);


	// Return true if string contains only spaces (as isspace() from <cctype>) or has zero length

	inline bool strempty(const char *str)

	{
		const char *ps=str;

		for( ; std::isspace(*ps); ps++)
			;
		return *ps == '\0';
	}

	inline bool strempty(const std::string& str)

	{
		return std::find_if(str.begin(),str.end(),std::not1(std::ptr_fun(is_space))) == str.end();
	}


	// Convert string to upper register 'in-place', returns pointer to itself

	inline char* strupper(char *str)

	{
		for(char *ps=str;*ps!='\0';ps++)
			*ps=std::toupper(*ps);
		return str;
	}


	inline std::string strupper(const std::string& str)

	{

		std::string str_new(str);

		for(int i = 0; i < str_new.size(); i++)
			str_new[i] = std::toupper(str_new[i]);

		return str_new;

	}


	//Return in-place pointer to c rightmost symbols in str.
	//If c>length(str), return str. If c==0, return pointer to "\0"
	//at the end of the str.

	//const char* strright(const char *str,int c);

	//inline char* strright(char *str,int c) {
	//  return const_cast<char*>(strright(const_cast<const char*>(str),c));
	//}


	//return true if string 's' ends with suffix 'sfx'

	//inline int isstrsfx(const char* s,const char* sfx)
	//{
	//   int csfx=std::strlen(sfx);
	//   return std::strcmp(strright(s,csfx),sfx)==0;
	//}

	/*

	// Return true if the argument string is not "tainted" in
	// approximately the same sence as Perl uses the term.
	// The second argument of enum TaintedType selects for what
	// kind of potentially tainted input to perform the check.
	// The function is modelled after the:
	// UNTAINT Perl module:
	// 19 September 1997 Rich Brown and Mark O'Neil
	// richard.e.brown@dartmouth.edu and mark.a.oneil@dartmouth.edu
	// Dartmouth College Computing Services
	//
	// In other words, this function checks if shell meta characters are
	// present in the string and returns false if this is the case.
	// Thus, 'not_tainted()' user-supplied strings are safer to pass as arguments
	// to functions which invoke subshell, such as popen() from standard C 
	// library.

	enum TaintedType { taintedFile, taintedPath, taintedEmail, taintedNumber, taintedShell };

	bool not_tainted(const char* s,UntaintType tainted_type);

	*/

} // namespace PRODDL

#endif // AT_STRING_UTIL_H__
