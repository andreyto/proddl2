//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_ARGPARSE_H__
#define AT_ARGPARSE_H__

// Helpers for parsing program arguments (on top of Boost program_options library).

#include "PRODDL/Common/string_util.hpp"
#include "PRODDL/Common/options.hpp"
#include <boost/program_options.hpp>
#include <boost/program_options/errors.hpp>
#include "PRODDL/Common/debug.hpp"
#include <exception>

namespace PRODDL {

	/// default conversion of argument name (as parsed in boost::program_options::variables_map)
	/// to options name (mirrors the deafult conversion by python argparse)
	std::string name_arg_to_opt(const std::string& name_arg) {
		return replace_char(name_arg,'-','_');
	}

	/// reciprocal of @see the name_arg_to_opt operation plus prepends "--"
	std::string name_opt_to_arg(const std::string& name_opt) {
		return "--"+replace_char(name_opt,'_','-');
	}

	/// throw required_option if argument is not defined
	inline void require_arg(
		const boost::program_options::variables_map& vm,
		const std::string name
		) throw (boost::program_options::required_option)
	{
			if( ! vm.count(name) ) {
				AT_THROW(boost::program_options::required_option(name));
			}
	}

	/// set key in Options from key in program arguments parsed into program_options::variables_map
	/// Options ket name is generated with @see name_arg_to_opt
	/// If name_arg is not found, does not set anything in opt
	template<typename T> std::string set_option_from_arg(
		const boost::program_options::variables_map& vm,
		Options& opt,
		const std::string& name_arg,
		bool require
		) {
			std::string name_opt = name_arg_to_opt(name_arg);
			if(vm.count(name_arg)) {
				opt.set(name_opt,vm[name_arg].as<T>());
			}
			else if(require) {
				AT_THROW(boost::program_options::required_option(name_arg));
			}
			return name_opt;
	}

	/// same as @see but accepts a default value instead of 'require' flag
	template<typename T> std::string set_option_from_arg_def(
		const boost::program_options::variables_map& vm,
		Options& opt,
		const std::string& name_arg,
		const T& defval 
		) {
			std::string name_opt = set_option_from_arg<T>(vm,opt,name_arg);
			if(! opt.has_option(name_opt) ) {
				opt.set(name_opt,defval);
			}
			return name_opt;
	}

} // namespace PRODDL

#endif // AT_ARGPRASE_H__
