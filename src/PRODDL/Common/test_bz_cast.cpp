//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/to_string.hpp"

#include "PRODDL/Common/bz_cast.hpp"

#include <blitz/array.h>

#include "gtest/gtest.h"

template<class T_num> void print_view(const blitz::Array<T_num,2>& a) {
  ATDBGOUT << "array" << a.shape() << "=" << std::endl;
  ATDBGOUT << '[' << std::endl;
  for(int i = a.lbound(0); i <= a.ubound(0); i++) {
    ATDBGOUT << '[' << std::endl;
    for(int j = a.lbound(1); j <= a.ubound(1); j++) {
      ATDBGOUT << a(i,j) << '(' << i << ',' << j << ')' << '\t';
    }
    ATDBGOUT << std::endl << ']' << std::endl;
  }
  ATDBGOUT << std::endl << ']' << std::endl;
}

template<class T_num> void print_view(const blitz::Array<T_num,3>& a) {
  ATDBGOUT << "array" << a.shape() << "=" << std::endl;
  ATDBGOUT << '[' << std::endl;
  for(int i = a.lbound(0); i <= a.ubound(0); i++) {
    ATDBGOUT << '[' << std::endl;
    for(int j = a.lbound(1); j <= a.ubound(1); j++) {
      ATDBGOUT << "[ ";
      for(int k = a.lbound(2); k <= a.ubound(2); k++) {
	ATDBGOUT << a(i,j,k) << '(' << i << ',' << j << ',' << k << ')' << '\t';
      }
      ATDBGOUT << " ]" << std::endl;
    }
    ATDBGOUT << std::endl << ']' << std::endl;
  }
  ATDBGOUT << std::endl << ']' << std::endl;
}

TEST(BzCastTest, All) {

  using namespace blitz;

  Array<TinyVector<int,2>,2> a(2,3,FortranArray<2>());
  TinyVector<int,2> x_ini(1,2);
  a = x_ini;

  ATDBGOUT << "a= " << std::endl;
  print_view(a);
  a.dumpStructureInformation(ATDBGOUT);
  ATOUTENDL();

  Array<int,3> a_control(2,3,2,FortranArray<3>());
  ATDBGOUT << "Control array of rank+1: " << std::endl;
  a_control.dumpStructureInformation(ATDBGOUT);  
  ATOUTENDL();

  Array<int,3> a_view = viewWithUnfoldedComponent<int>(a,0,1);

  ATDBGOUT << "a_view= " << std::endl;
  print_view(a_view);
  a_view.dumpStructureInformation(ATDBGOUT);
  ATOUTENDL();

  Array<int,3> a_view2 = viewWithUnfoldedComponent<int>(a,1,4);
  // Comment out the line above and uncomment the line below: you must get compilation
  // error since sizeof(TinyVector<int,2>) % sizeof(TinyVector<char,3>) != 0
  //Array<TinyVector<char,3>,3> a_view2 = viewWithUnfoldedComponent<TinyVector<char,3> >(a,1,4);

  ATDBGOUT << "a_view2= " << std::endl;
  print_view(a_view2);
  a_view2.dumpStructureInformation(ATDBGOUT);
  ATOUTENDL();

  Array<TinyVector<int,2>,2> a_view_orig = viewWithFoldedComponent<TinyVector<int,2> >(a_view2);
  ATDBGOUT << "a_view_orig= " << std::endl;
  print_view(a_view_orig);
  a_view_orig.dumpStructureInformation(ATDBGOUT);
  ATOUTENDL();
}

