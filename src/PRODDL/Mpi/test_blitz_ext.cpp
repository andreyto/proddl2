//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//


#include <iostream>

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Mpi/blitz.hpp"

int
main_test()
{
    using namespace std;
    using namespace PRODDL::Mpi;
  cout << "main_test(): Going to call OOMPI_COMM_WORLD.Init()" << endl;
  int argc = 0;
  char **argv = 0;
  OOMPI_COMM_WORLD.Init(argc,argv);
  
  cout << "main_test(): OOMPI_COMM_WORLD.Init() called" << endl;
  
  int rank = OOMPI_COMM_WORLD.Rank();
  int size = OOMPI_COMM_WORLD.Size();
  
  cout << "Blitz MPI datatypes testing. I am " << rank << " of " << size << endl;

  typedef blitz::TinyVector<double,3> Point;
  Point x;

  typedef blitz::Array<double,2> ArrayNum;
  ArrayNum y(2,2);

  typedef blitz::Array<Point,2> ArrayPoint;
  ArrayPoint z(2,2);

  typedef ArrayMessenger<Point,2> ArrMsgr;

  if(rank == 0) {
    x = 1,2,3;
    y = 0,1,2,3;
    z = Point(1,1,1), Point(2,2,2), Point(3,3,3), Point(4,4,4);
    for(int i = 1; i < size; i++) {
      OOMPI_COMM_WORLD[i].Send(message(x));
      OOMPI_COMM_WORLD[i].Send(message(y));
      OOMPI_COMM_WORLD[i].Send(message(z));
      cout << "Process with rank " << rank << " sent data " << x << "  " << y << "  " << z <<endl;
      z = Point(i);
      ArrMsgr arr_msgr(z);
      OOMPI_COMM_WORLD[i].Send(arr_msgr.messageStruct());
      OOMPI_COMM_WORLD[i].Send(arr_msgr.messageData());
      cout << "Process with rank " << rank << " sent data as ArrayMessenger(z)" << z << endl;
      ArrayPoint z_sub(z(blitz::Range(0,1),blitz::Range(0,0)));
      ArrMsgr arr_msgr_sub(z_sub);
      OOMPI_COMM_WORLD[i].Send(arr_msgr_sub.messageData());
      cout << "Process with rank " << rank << " sent data as ArrayMessenger(z_sub)" << z_sub << endl;
    }
  }
  else {
    OOMPI_COMM_WORLD[0].Recv(message(x));
    OOMPI_COMM_WORLD[0].Recv(message(y));
    OOMPI_COMM_WORLD[0].Recv(message(z));
    cout << "Process with rank " << rank << " received data " << x << " " << y << " " << z << endl;
    ArrayPoint z0;
    ArrMsgr arr_msgr(z0);
    OOMPI_COMM_WORLD[0].Recv(arr_msgr.messageStruct());
    OOMPI_Status status = OOMPI_COMM_WORLD[0].Recv(arr_msgr.messageData());
    cout << "Process with rank " << rank << " received data as ArrayMessenger(z0)" << z0 << endl;
    OOMPI_Message msg = arr_msgr.messageData();
    ATOUTVAR(msg.Get_count()); ATOUTVAR(z0.size()); ATOUTVAR(sizeof(z0(0,0))); ATOUTENDL();
    ATOUTVAR(status.Get_count(OOMPI_Datatype(msg.Get_type()))); ATOUTENDL();
    status = OOMPI_COMM_WORLD[0].Recv(arr_msgr.messageData());
    cout << "Process with rank " << rank << " received data as ArrayMessenger(z0)" << z0 << endl;
    ATOUTVAR(msg.Get_count()); ATOUTVAR(z0.size()); ATOUTVAR(sizeof(z0(0,0))); ATOUTENDL();
    ATOUTVAR(status.Get_count(OOMPI_Datatype(msg.Get_type()))); ATOUTENDL();
    ATOUTVAR(elements(status,msg,z0));  ATOUTENDL();
  }

  OOMPI_COMM_WORLD.Finalize();

  return 0;
}
  
#include "PRODDL/Testing/ctestc.hpp"

CTEST_MODULE()
