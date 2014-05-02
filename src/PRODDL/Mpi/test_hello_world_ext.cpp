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
  cout << "main_test(): Going to call OOMPI_COMM_WORLD.Init()" << endl;
  int argc = 0;
  char **argv = 0;
  OOMPI_COMM_WORLD.Init(argc,argv);
  
  cout << "main_test(): OOMPI_COMM_WORLD.Init() called" << endl;
  
  int rank = OOMPI_COMM_WORLD.Rank();
  int size = OOMPI_COMM_WORLD.Size();
  
  cout << "Hello World! I am " << rank << " of " << size << endl;

  OOMPI_COMM_WORLD.Finalize();

  return 0;
}
  
#include "PRODDL/Testing/ctestc.hpp"

CTEST_MODULE()
