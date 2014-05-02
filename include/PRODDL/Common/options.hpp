//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_OPTIONS_H__
#define PRODDL_OPTIONS_H__

/*
Class 'Options' to represent optional parameters.
It has an interface of a polymorphic dictionary.
The typechecking is performed during extraction.
Each 'Options' object can in turn contain other
named 'Options' objects ("blocks of options").
Example:
...
Options options;
int main() {
options.set("epsilon",1e-7);
options.set("max_iter",1000);
Options ff_options;
Options &ff_options = options.setBlock("ForceField");
ff_options.set("lj_eps",0.3);
ff_options.set("lj_r",0.37");
...
}
...
void optimize() {
extern Options options;
float epsilon;
options.get("epsilon",epsilon);
float lj_eps;
options.getBlock("ForceField").get("lj_eps",lj_eps);
...
}
...

This class has a deep copy semantics.
Requirements on the value types of options:
They must be "simple" types: copy constructable, copyable, constructor
should not throw exceptions.
Implementation uses C++ typeid.
*/

#include <boost/any.hpp>

#include <string>

#include <map>

#include "PRODDL/exceptions.hpp"

#include <iostream>

namespace PRODDL {


	class options_key_error : public docktk_error {
	public:
		options_key_error(const std::string& msg)  throw (): 
		  docktk_error(msg) {
		  }
		  virtual ~options_key_error() throw (){}
	};



	class Options {

	public:

		typedef std::string KeyType;

		typedef Options Self;

	public:

		template<typename ValueType> 
		void
			get(const KeyType& key, ValueType& value) const {

				MapLeavesType::const_iterator iter = map_leaves.find(key);
				if( iter == map_leaves.end() )
					throw options_key_error(key);
				value = boost::any_cast<ValueType>(iter->second);

		}

		template<typename ValueType, typename DefValueType> 
		void
			getdefault(const KeyType& key, ValueType& value, const DefValueType& default_value) const {
				if( has_option(key) )
					get(key,value);
				else
					value = default_value;
		}


		bool has_option(const KeyType& key) const {
			return ! ( map_leaves.find(key) == map_leaves.end() );
		}

		bool has_block(const KeyType& key) const {
			return ! ( map_blocks.find(key) == map_blocks.end() );
		}

		template<typename ValueType>
		ValueType&
			set(const KeyType& key, const ValueType& value) {
				map_leaves[key] = boost::any(value);
				return *boost::any_cast<ValueType>(&map_leaves[key]);
		}

		Self&
			set(const KeyType& key, const Self& value) {
				map_blocks[key] = value;
				return map_blocks[key];
		}

		void
			remove(const KeyType& key) {
				map_leaves.erase(key);
				map_blocks.erase(key);
		}

		Self&
			getBlock(const KeyType& key) {
				MapBlocksType::iterator iter = map_blocks.find(key);
				if( iter == map_blocks.end() )
					throw options_key_error(key);
				return iter->second;
		}

		const Self&
			getBlock(const KeyType& key) const {
				MapBlocksType::const_iterator iter = map_blocks.find(key);
				if( iter == map_blocks.end() )
					throw options_key_error(key);
				return iter->second;
		}


		void clear() {
			map_leaves.clear();
			map_blocks.clear();
		}

		// debugging output

		void
			print(std::ostream& out) const {

				for(MapLeavesType::const_iterator p = map_leaves.begin();
					p != map_leaves.end();
					p++) {
						out << p->first << "\n";
				}
				out << "Options object.";

		}

	protected:

		typedef std::map<KeyType,boost::any> MapLeavesType;

		typedef std::map<KeyType,Self> MapBlocksType;



		MapLeavesType map_leaves;

		MapBlocksType map_blocks;

	};

	inline
		std::ostream&
		operator<< (std::ostream& out, const Options& opts) {

			opts.print(out);
			return out;

	}


} // namespace PRODDL

#endif // PRODDL_OPTIONS_H__
