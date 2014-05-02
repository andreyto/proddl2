#include "PRODDL/IO/hdf5.hpp"
#include "PRODDL/Common/bz_cast.hpp"
#include "gtest/gtest.h"
#include <boost/filesystem.hpp>

// Fixture
class BobHdf5Test : public ::testing::Test {
protected:
	std::string filename;
	PRODDL::Hdf5::HDF5File *p_file;
public:
	virtual void SetUp() {
		filename = "tmp_test_bob_hdf5.hdf5";
		p_file = 0;
		p_file = new PRODDL::Hdf5::HDF5File(filename,PRODDL::Hdf5::HDF5File::inout);
	}

	virtual void TearDown() {
		delete p_file;
		boost::filesystem::remove(filename);
	}

};

TEST_F(BobHdf5Test, Array2D) {
	blitz::Array<double,2> a(4,2);
	a = 1, 2, 
		3, 4,
		5, 6,
		7, 8;

	blitz::Array<double,2> at = a.transpose(1,0);
	p_file->setArray("at", at);

	// Read it and compare to original
	blitz::Array<double,2> at_read;
	at_read.reference(p_file->readArray<double,2>("at"));
	EXPECT_TRUE(blitz::all(at==at_read));

	// Read part of it
	//blitz::Array<double,2> at_read1(2,2);
	blitz::Array<double,1> at_read1(4);
	p_file->readArray("at",0,at_read1);
	EXPECT_TRUE(at_read1.size() == 4);
	EXPECT_TRUE(blitz::all(at(0,blitz::Range::all())==at_read1));

	blitz::Array<double,1> at_read2(4);
	at_read2 = 0;

	for(int i=0; i<4; ++i) {
		p_file->appendArray("at2",at_read1);
	}

	p_file->readArray("at2",2,at_read2);
	EXPECT_TRUE(at_read2.size() == 4);
	EXPECT_TRUE(blitz::all(at_read1==at_read2));

}

TEST_F(BobHdf5Test, Array1DPoints) {
	using namespace blitz;
	typedef TinyVector<double,2> Point;
	Point p1(1,1);
	Point p2(2,2);
	Array<Point,1> a(1000000);
	a(0) = p1;
	a(1) = p2;
	Array<double,2> a_f(viewWithFlattenedComponent(a));
	for(int i=0; i<1; ++i) {
		p_file->appendArray("a_f", a_f);
	}
}
