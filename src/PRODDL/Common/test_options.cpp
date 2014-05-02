//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Common/options.hpp"
#include "PRODDL/Common/options_io_json.hpp"

#include "PRODDL/Common/debug.hpp"

#include "gtest/gtest.h"

using namespace PRODDL;

namespace {

	Options opts;
	
	void test_set() {
		opts.set("epsilon",1e-3);
		opts.set("niter",1000);
		Options& block = opts.set("lj",Options());

		block.set("r",0.75);
		block.set("e",0.34);
	}

	void test_get() {

		double epsilon;
		opts.get("epsilon",epsilon);
		ATOUTVAR(epsilon);
		ATALWAYS(epsilon==1e-3,"Unexpected option's value");

		int niter;
		opts.getdefault("niter",niter,10);
		ATOUTVAR(niter);
		ATALWAYS(niter == 1000,"Unexpected option's value");

		Options& block = opts.getBlock("lj");
		double r;
		block.get("r",r);
		ATOUTVAR(r);
		ATALWAYS(r == 0.75,"Unexpected option's value");

		block.remove("e");

		try {

			double e;
			block.get("e",e);

		}
		catch(options_key_error exception_msg) {
			ATOUTVAR(exception_msg.what()); ATOUTENDL();
		}

	}
	
	const char inp_json_str[] = 
		"{\"epsilon\":1e-3,"
		" \"niter\":1000,"
		" \"lj\":{\"r\":0.75,\"e\":0.34}"
		"}";

	void test_json_io() {
		opts.clear();
		test_set();
		load_options_from_json_string(inp_json_str,opts);
		test_get();
	}

}

TEST(OptionsTest, All) {

	test_set();
	test_get();
	test_json_io();
}

