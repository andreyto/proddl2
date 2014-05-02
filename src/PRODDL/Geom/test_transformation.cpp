//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/common_algor.hpp"

#include "PRODDL/Geom/transformation.hpp"

#include "PRODDL/Testing/exception.hpp"

#include "gtest/gtest.h"

TEST(TransformationTest, All) {

  using namespace PRODDL::Geom;
  using namespace PRODDL::Math;
  typedef double T_num;
  typedef Translation<T_num> Tran;
  typedef Rotation<T_num> Rot;
  typedef RotationTranslation<T_num> RotTran;
  typedef Tran::Point Point;
  typedef Tran::Matrix3 Matr3;
  Point cm, v_trans;
  cm = 5,4,3;
  v_trans = 1,2,3;
  Matr3 matr_rot, matr_rot2;
  matr_rot = 
    0,1,1,
    0.2,0.3,0.4,
    1,0.3,0.5;
  matr_rot2 = 
    0.4,1.7,1.8,
    0.7,0.9,0.3,
    1.3,0.1,0.9;
  RotTran rt =  Rot(matr_rot2)*Tran(cm + v_trans)*(Rot(matr_rot)*Tran(0-cm));
  rt = rt.inverse();
  ATOUTVAR(rt.getTensor()); ATOUTVAR(rt.getVector()); ATOUTENDL();
  Matr3 matr_rot_res;
  matr_rot_res = 
    2.14,0.48,0.92,
    1.45,1.06,1.6,
    1.98,1.21,1.79;
  Point vect_res(-0.2324,1.5222,1.3322);
  if( ! ( are_all_close(rt.getTensor(),matr_rot_res) && are_all_close(rt.getVector(),vect_res)) ) {
    throw PRODDL::Testing::test_error("Transform: Test output does not match the expected values.");
  }

  // test operator()()
  Matr3 matr_rot_rot;
  matr_rot_rot = rt.rotation()(matr_rot);
  ATOUTVAR(matr_rot_rot); ATOUTENDL();
  Matr3 matr_rot_rot_res;
  matr_rot_rot_res = 
    3.61830443,   1.20609827,   2.05409441,
    17.35934489,  -1.37593931,   0.66919557,
    -13.74527938,   0.37978324,  -1.44236513;
  if( ! are_all_close(matr_rot_rot,matr_rot_rot_res) ) {
    throw PRODDL::Testing::test_error("Transform: Test output does not match the expected values.");
  }
}

