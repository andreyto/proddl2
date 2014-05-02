//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Common/string_util.hpp"


namespace PRODDL {

std::string strtranc(const std::string& s)
{
   std::string::const_iterator pb = std::find_if(s.begin(),s.end(),std::not1(std::ptr_fun(is_space)));
   std::string::const_iterator pe = 
      std::find_if(s.rbegin(),
      std::string::const_reverse_iterator(pb),
      std::not1(std::ptr_fun(is_space))).base();
   std::string sret(pb,pe);
   return sret;
}

  /*
char *strtranc(char *dest, const char *src)

{

   const char *psbeg=src;
   for(;std::isspace(*psbeg);psbeg++)

   ;

   if(*psbeg)

   {

      const char *psend=src+std::strlen(src);

      while(std::isspace(*--psend));

      int len;

      std::strncpy(dest,psbeg,len=psend-psbeg+1)[len]=0;

   }

   else

   *dest='\0';

   return dest;

}



const char* strright(const char *str,int c)
{
   const char *p=str;
   for(;*p!='\0';p++)
      ;
   const char *pleft=p-c;
   if(pleft<str)
      pleft=str;
   return (char*) pleft;
}

  */

/*

bool not_tainted(const char* s,UntaintType tainted_type) {
  char * reg_str = 0;
  bool must_match;
  if (tainted_type == taintedFile) {
    reg_str = "[\\w\\-._]+";  // includes a-zA-Z0-9._ and -
    must_match = true;
  }
  else if(tainted_type == taintedPath) {
    reg_str = "[\\/\\w\\-._]+"; // includes 'file' chars + "/"
    must_match = true;
  }
  // TODO: e-mail requires less restrictive pattern than inplemented,
  // see Perl::Mail::CheckUser.pm
  else if(tainted_type == taintedEmail) {
    reg_str = "[\\w\\-'._]+\\@[\\w\\.-]+"; // a-zA-Z0-9.'_- then @ then a domain name
    must_match = true;
  }
  else if(tainted_type ==  taintedNumber) {
    reg_str = "[0-9]*\\.?[0-9]*"; //  ddd.ddd
    must_match = true;
  }
  else if(tainted_type == taintedShell) {
    reg_str = "[\\&\\;\\`\\'\\\\\"\\|\\*\\?\\~\\<\\>\\^\\(\\)\\[\\]\\{\\}\\$\\n\\r]"; 
    must_match = false;
  }
  // Now create regular expression and try to match
}

*/

} // namespace PRODDL
