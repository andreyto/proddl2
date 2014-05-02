//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/docking.hpp"
#include "gtest/gtest.h"

using namespace PRODDL;
using namespace std;

typedef double T_num;
typedef Docking<T_num> D;

class DockIO_BinTest : public ::testing::Test {

protected:

	D::Rotation rot; //identity

	D::TranValues tv;

public:

	virtual void SetUp() {

		tv.reference(D::TranValues(2));
		tv(0).value = 1.;
		tv(0).tran = D::Translation(1.);

		tv(1).value = 2.;
		tv(0).tran = D::Translation(2.);

	}
};

TEST_F(DockIO_BinTest, RotFftScanIO_Bin) {

	{
		D::RotFftScanIO_Bin io("dock_rot_scan.tmp.bin",ios::out|ios::trunc);
		io.write_record(rot,tv);
		tv(1).value = 3.;
		io.write_record(rot,tv);
	}
	{
		D::TranValues tv_in;
		D::RotFftScanIO_Bin io("dock_rot_scan.tmp.bin",ios::in);
		bool status=false;
		status = io.read_record(rot,tv_in);
		EXPECT_TRUE(status);
		EXPECT_TRUE(tv_in.size() == 2 && tv_in(1).value==2.);
		status = io.read_record(rot,tv_in);
		EXPECT_TRUE(status);
		EXPECT_TRUE(tv_in.size() == 2 && tv_in(1).value==3.);
		status = io.read_record(rot,tv_in);
		EXPECT_FALSE(status);
		status = io.read_record(rot,tv_in);
		EXPECT_FALSE(status);

	}
}

TEST_F(DockIO_BinTest, RotFftScanIO_Collector) {

	{
		D::RotFftScanIO_Bin io("dock_rot_scan.tmp.1.bin",ios::out);
		io.write_record(rot,tv);
		tv(1).value = 3.;
		io.write_record(rot,tv);
	}
	{
		D::RotFftScanIO_Bin io("dock_rot_scan.tmp.2.bin",ios::out);
		io.write_record(rot,tv);
		tv(1).value = 4.;
		io.write_record(rot,tv);
	}
	{
		ofstream task_list("dock_rot_scan.tasks.tmp.tab");
		task_list << "dock_rot_scan.tmp.1.bin" << "\n";
		task_list << "dock_rot_scan.tmp.2.bin" << "\n";
	}
	{
		D::TranValues tv_in;
		D::RotFftScanIO_Collector<D::RotFftScanIO_Bin> task_coll("dock_rot_scan.tasks.tmp.tab");
		bool status=false;
		for(int i=0; i<4; ++i) {
			status = task_coll.read_record(rot,tv_in);
			EXPECT_TRUE(status);
			if(i==1 || i==2)
				EXPECT_TRUE(tv_in.size() == 2 && tv_in(1).value==3.);
			if(i==3)
				EXPECT_TRUE(tv_in.size() == 2 && tv_in(1).value==4.);
		}
		status = task_coll.read_record(rot,tv_in);
		EXPECT_FALSE(status);
		status = task_coll.read_record(rot,tv_in);
		EXPECT_FALSE(status);

	}
}
