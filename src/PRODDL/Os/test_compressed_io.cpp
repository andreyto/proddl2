//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Os/compressed_io.hpp"

#include <string>

#include "PRODDL/Testing/exception.hpp"

#include "PRODDL/Testing/config.hpp"

#include "PRODDL/Common/to_string.hpp"

#include "gtest/gtest.h"

TEST(CompressedIoTest, All) {  
	const char* extensions[] = {".Z", ".gz"};
	for(int i = 0; i < sizeof(extensions)/sizeof(extensions[0]); i++) {
		std::string fileName = std::string("test_compressed_io")+extensions[i];
		boost::shared_ptr<std::ostream> p_out = PRODDL::openCompressedOutput(fileName.c_str());
		int j;
		for(j = 0; j < 100; j++)
			(*p_out) << j << std::endl;
		(*p_out) << j;
		p_out.reset();
		boost::shared_ptr<std::istream> p_in = PRODDL::openUncompressedInput(fileName.c_str());
		for(j = 0; j <= 100; j++) {
			int j_received = -1;
			(*p_in) >> j_received;
			std::cout << "expected: " << j << "   received: " << j_received << std::endl;
			if( j != j_received ) {
				throw PRODDL::Testing::test_error(std::string("input does not match what was written during output. Expected: ") + 
					PRODDL::to_string(j) + ". Received: " + PRODDL::to_string(j_received));
			}
		}
		if( ! p_in->eof() ) {
			throw PRODDL::Testing::test_error("Expected eof() in input uncompressed stream");
		}
	}
}
