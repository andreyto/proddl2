//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_REGEXP_H__
#define PRODDL_REGEXP_H__

#ifdef GNURX
   #include <rxposix.h> //GNU Posix regular expressions from GNU rx lib
#else
   #include <regex.h> //Posix regular expressions supplied by system
#endif //ifdef GNURX

#include <string>
#include <vector>
#include <algorithm>


//wrapper class around POSIX regular expressions library.
//Example of usage:
/*
   //Parse assignment of kind 'y = x,z' and print found arguments and entire expression.
   RegExp rx("pattern[[:blank:]]*=[[:blank:]]*\\([[:alnum:]]"
   "\\{1,\\}\\)[[:blank:]]*,[[:blank:]]*\\([[:alnum:]]*\\)",REG_NEWLINE);
   cout << "Looking for pattern: " << rx.GetPatt() << endl;	
   if(rx.good())
	{
  	  RegExp::astr avalues;
     string text("pattern = MyPattern, 2\n Something Else");
	  cout << "Text is: " << text << endl;
      if(rx.split(text,avalues) == 0)
      {
         for(size_t i = 1; i < avalues.size(); i++)
         {
            cout << "Subexpression #" << i << " is: " << avalues[i] << endl;
         }
         string sentire = rx.extract(text,0);
         cout << "Entire match is: " 
            << sentire << endl;
      }
      else
      {
         cout << "Nothing found." << endl;
      }
   }
   else
   {
      cout << "Regular expression compilation error: " 
         << rx.GetErrorMsg() << endl;
   }
*/

namespace PRODDL {

class RegExp
{
   public:
   
   typedef std::vector<std::string> astr;
   typedef regmatch_t mat;
   
   protected:
   
   std::string spatt; //pattern string
   int cflags; //compilation flags
   regex_t r_t; //compiled reg. expr.
   //array with offsets of found subexpressions. Number
   //of elements is r_t.re_nsub.
   regmatch_t *ar_m;
   int car_m; //size of ar_m
   int state; //error code returned by regcomp()
   std::string errmsg; //error message obtained from regerror()
   
   void get_regerror()
   {
      size_t length = regerror(state,&r_t,NULL,0);
      char *buffer = new char[length];
      buffer[0]='\0';
      regerror(state,&r_t,buffer,length);
      errmsg = buffer;
      delete[] buffer;
   }

   //Compile pattern spatt to reg. expr. r_t.
   //Resources (r_t and ar_m) must be freed by freeres() 
   //or never been allocated 
   //before calling this function.
   void compile()
   {
      state = regcomp(&r_t,spatt.c_str(),cflags);
      if(state != 0) //compilation error
      {
         get_regerror();
         car_m = 0;
         ar_m = 0;
      }
      else //compiled OK
      {
         errmsg = "";
         car_m = r_t.re_nsub+1;
         ar_m = new regmatch_t[car_m];
      }
   }

   //Frees resources. Compile must be called before
   //calling this function.
   void freeres()
   {
      regfree(&r_t);
      delete[] ar_m;
   }
   
   //service function used by overloaded versions of split().
   //See split() description.
   void _split_t(const char *_str,RegExp::astr& asubs) const
   {
      asubs.clear();
      asubs.resize(car_m);
      for(size_t i = 0; i != asubs.size(); i++)
      {
         const regmatch_t& m = ar_m[i];
         if(m.rm_so != -1 && m.rm_eo != -1)
         {
            asubs[i].append(_str+m.rm_so,_str+m.rm_eo);
         }
      }
   }

   public:
   
   explicit RegExp(const std::string& _spatt,int _cflags=REG_NEWLINE|REG_EXTENDED) :
      spatt(_spatt), cflags(_cflags)
   {
      compile();
   }
   explicit RegExp(const char* _spatt,int _cflags=REG_NEWLINE|REG_EXTENDED) :
      spatt(_spatt), cflags(_cflags)
   {
      compile();
   }
   
   RegExp(const RegExp& x) :
      spatt(x.spatt), cflags(x.cflags)
   {
      compile();
   }
   RegExp& operator= (const RegExp& x)
   {
      if( this != &x )
      {
         freeres();
         spatt = x.spatt;
         cflags = x.cflags;
         compile();
         //paranoid safety feature:
         int c = car_m < x.car_m ? car_m : x.car_m;
         std::copy(x.ar_m,x.ar_m+c,ar_m);
      }
      return *this;
   }
   ~RegExp()
   {
      freeres();
   }
   //search for reg. expr. in string _str.
   //From regex man page:
   //Non zero return values mean; 
   //REG_NOMATCH - not found, REG_ESPACE - out of memory.
   //_eflags can be: REG_NOTBOL, REG_NOTEOL.
   //We provide overloaded versions for 
   //{const char*|const string&}.
   int match(const char *_str,int _eflags=0)
   {
      return regexec(&r_t,const_cast<char*>(_str),car_m,ar_m,_eflags);
   }
   int match(const std::string& _str,int _eflags=0)
   {
      return match(_str.c_str());
   }
   //call match and fill array of strings with subexpressions pointed by
   //ar_m entries as set by regexec(). 
   //Return result of match(). If the entire expression is not found, this
   //function does not change the contents of asubs argument.
   //Remember, that in ar_m the first element (with offset 0) 
   //points to the entire matched regular expression. 
   //We use the same scheme for 'asubs' argument.
   //Note, that this function does not distinguish between the alternatives when
   //subexpression matched empty string or was not used at all. In both
   //cases the result will be empty string in corresponding element of asubs.
   //We provide overloaded versions for
   //{const char*|const string&}.
   int split(const char *_str,RegExp::astr& asubs,int _eflags=0)
   {
      int ret = match(_str,_eflags);
      if( ret == 0) _split_t(_str,asubs);
      return ret;
   }
   int split(const std::string& _str,RegExp::astr& asubs,int _eflags=0)
   {
      const char *__s = _str.c_str();
      int ret = match(__s,_eflags);
      if( ret == 0) _split_t(__s,asubs);
      return ret;
   }
   const std::string& GetErrorMsg() { return errmsg;  } const
   int GetState() { return state; } const
   int ClearState()
   {
      int _t = state;
      state = 0;
      errmsg = "";
      return _t;
   }
   bool good() const { return state==0; }
   const std::string& GetPatt() const { return spatt; }
   //return array of matches (1-st element corresponds to the entire
   //expression:
   const RegExp::mat* GetMatches() const { return ar_m; } 
   //return i-th element of array of matches:
   const RegExp::mat& GetMatch(int i) const { return ar_m[i]; }
   //return string pointed by i-th element of ar_m, i.e.
   //i-th subexpression. i=0 corresponds to entire expression.
   std::string extract(const char *_str,int i) const
   {
      const regmatch_t& m = ar_m[i];
      if(m.rm_so != -1 && m.rm_eo != -1)
      {
         return std::string(_str+m.rm_so,_str+m.rm_eo);
      }
      return std::string();
   }
   std::string extract(const std::string& _str, int i) const
   {
      return extract(_str.c_str(),i);
   }
   //Return ponter to the beginning of i-th
   //match. _str argument is supposed to be the same as
   //used in last call to match(). If i-th subexpression
   //did not match, return 0.
   const char* BegMatch(const char *_str,int i) const
   {
      const regmatch_t& m = ar_m[i];
      if(m.rm_so != -1)
         return _str+m.rm_so;
      return 0;
   }
   //Return ponter to the "past the end" position of i-th
   //match. _str argument is supposed to be the same as
   //used in last call to match(). If i-th subexpression
   //did not match, return 0.
   const char* EndMatch(const char *_str,int i) const
   {
      const regmatch_t& m = ar_m[i];
      if(m.rm_eo != -1)
         return _str+m.rm_eo;
      return 0;
   }
   //return number of parenthised subexpressions plus one 
   //(plus one is because first element
   //in ar_m is used for the match of entire expression):
   int CntSub() const { return car_m; }
};


inline std::string RxGetParam(const char *spatt,const char *text,int isub=1)
{
   RegExp rx(spatt,REG_EXTENDED);
   if(rx.match(text)==0 && rx.CntSub()>isub)
   {
      return rx.extract(text,isub);
   }
   return std::string();
}

inline const char* RxFindBeg(const char *spatt,const char *text,int isub=0)
{
   RegExp rx(spatt,REG_EXTENDED);
   if(rx.match(text)==0 && rx.CntSub()>isub)
   {
      return rx.BegMatch(text,isub);
   }
   return 0;
}

inline const char* RxFindEnd(const char *spatt,const char *text,int isub=0)
{
   RegExp rx(spatt,REG_EXTENDED);
   if(rx.match(text)==0 && rx.CntSub()>isub)
   {
      return rx.EndMatch(text,isub);
   }
   return 0;
}


} // namespace PRODDL

#endif // DOKTK_REGEXP_H__
