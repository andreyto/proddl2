//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_OPTIONS_IO_JSON_H__
#define PRODDL_OPTIONS_IO_JSON_H__

/*
*
Methods to load Options class from JSON string or file.
*
*/

#include "PRODDL/Common/options.hpp"

namespace PRODDL {

	void load_options_from_json_file(std::string inp_file,Options& opt);
	
	void load_options_from_json_string(std::string inp_str,Options& opt);

} // namespace PRODDL

#endif // PRODDL_OPTIONS_IO_JSON_H__
