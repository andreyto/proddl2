//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_MPI_DISPATCHER_H__
#define PRODDL_MPI_DISPATCHER_H__

#include "PRODDL/Mpi/mpi.hpp"

#include "PRODDL/Mpi/tags.hpp"

#include <unordered_map>

namespace PRODDL {

  namespace Mpi {

template < class HandlerClass >
class Dispatcher {

public:

  typedef void (HandlerClass::* HandlerMethod) (OOMPI_Status&);

protected:

  typedef unordered_map<int, HandlerMethod> HandlerMap;

public:

  void setHandlerObject(HandlerClass& x) {

    p_handler = &x;

  }

  void setHandlerMethod(int tag, HandlerMethod method) {

    handlerMap[tag] = method;

  }

  void removeHandlerMethod(int tag) {

    handlerMap.erase(tag);

  }


  virtual
  void run() = 0;

protected:

  void handleMessage(int tag, OOMPI_Status& status) {

    

  }

protected:

  HandlerClass *p_handler;

  HandlerMap handlerMap;


}; // class Dispatcher


template < class HandlerClass >
class SlaveDispatcher : public Dispatcher<HandlerClass> {

public:

  typedef Dispatcher<HandlerClass> Base;

  virtual
  void run() {

    this->handleMessage(tag,status);

  }

}; // class SlaveDispatcher


template < class HandlerClass >
class MasterDispatcher : public Dispatcher<HandlerClass> {

public:

  typedef Dispatcher<HandlerClass> Base;


}; // class MasterDispatcher


  }} // namespace PRODDL::Mpi

#endif // PRODDL_MPI_DISPATCHER_H__
