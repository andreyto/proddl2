//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//


#include <iostream>
#include "oompi.h"

int
main_test()
{

    using namespace std;

  // Init() and Finalize() are moved to a separate module (and oompi was compiled as a shared library)
  //cout << "main_test(): Going to call OOMPI_COMM_WORLD.Init()" << endl;
  //OOMPI_COMM_WORLD.Init();
  
  //cout << "main_test(): OOMPI_COMM_WORLD.Init() called" << endl;
  
  int rank = OOMPI_COMM_WORLD.Rank();
  int size = OOMPI_COMM_WORLD.Size();
  
  cout << "Hello World! I am " << rank << " of " << size << endl;

  //OOMPI_COMM_WORLD.Finalize();

  return 0;
}
  
#include "PRODDL/Testing/ctestc.hpp"

CTEST_MODULE()
