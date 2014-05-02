//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_COMMON_QUEUE_H__
#define PRODDL_COMMON_QUEUE_H__

// Defines BoundPriorityQueue - a priority queue that has a restriction on its size
// (the top element will be popped away once upon a new element insertion once a size
// limit is reached) and a restriction on element value (elements larger than a certain
// value will never be inserted).
// Implemented similar to std::priority_queue, with additional constraints desribed above,
// and also with access to the internal container.
// Access to internal container is provided for the benefit of bulk operations such as
// MPI sends.
// Requirements to the internal container type: those listed for std::priority_queue,
// plus a 'reserve' function (for performance reasons). std::vector<> is a default container type.



#include <functional>

#include <vector>

#include <limits>

#include "PRODDL/Common/logger.hpp"

#include "PRODDL/exceptions.hpp"

namespace PRODDL {


  class BoundPriorityQueueError : public docktk_error {
  public:
    BoundPriorityQueueError(const std::string& msg)  throw (): 
      docktk_error(msg) {
    }
    virtual ~BoundPriorityQueueError() throw (){}
  };


  ///////////////////////////////////////////////////////////

  template <class _Tp, 
	    class _Sequence = std::vector<_Tp>,
	    class _Compare  = std::less<typename _Sequence::value_type>,
	    class _LimitType = _Tp,
	    class _LimitCompare = _Compare>

  class BoundPriorityQueue
  {

    typedef typename _Sequence::value_type _Sequence_value_type;


  public:
    typedef typename _Sequence::value_type      value_type;
    typedef typename _Sequence::size_type       size_type;
    typedef          _Sequence                  container_type;

    typedef          _LimitType                 limit_type;

    typedef typename _Sequence::reference       reference;
    typedef typename _Sequence::const_reference const_reference;


  protected:

    _Sequence c;
    _Compare comp;
    _LimitCompare limit_comp;

    size_type size_limit;

    limit_type max_val_limit;

    bool is_valid;

  public:

    BoundPriorityQueue() :
      size_limit(0),
      is_valid(false)
    {}

    explicit BoundPriorityQueue(size_type _size_limit,
				const limit_type& _max_val_limit,
				const _Compare& __x = _Compare(),
				const _LimitCompare& __y = _LimitCompare()) {
      init(_size_limit,_max_val_limit,__x,__y);
    }

    void init(size_type _size_limit,
	      const limit_type& _max_val_limit,
	      const _Compare& __x = _Compare(),
	      const _LimitCompare& __y = _LimitCompare()) {

      size_limit = _size_limit;
      //Why would I write this?:
      //if(size_limit == 0) size_limit = std::numeric_limits<size_type>::max();    
      ATLOG_ASSERT_1(size_limit > 0);
      max_val_limit = _max_val_limit;
      comp = __x;
      limit_comp = __y;
      clear();
    }


    bool empty() const { return c.empty(); }

    size_type size() const { return c.size(); }

    const_reference top() const { 

      ATLOG_ASSERT_1( is_valid );

      //if( ! is_valid )
      //throw BoundPriorityQueueError("Invalid state for top/push/pop operations");

      return c.front(); 

    }

    void clear() {

      ATLOG_TRACE_4;

      c.clear();
      if( size_limit ) {
	c.reserve(size_limit);
	is_valid = true;
      }

    }


    const _Sequence& data() const {

      return c;

    }

    _Sequence& data() {

      return c;

    }


    const _Compare& getCompare() const {

      return comp;

    }

    _Compare& getCompare() {

      return comp;

    }

    const _LimitCompare& getLimitCompare() const {

      return limit_comp;

    }

    _LimitCompare& getLimitCompare() {

      return limit_comp;

    }

    // Top, push and pop cannot be called any more after the call to this function
    // (because internal sequence is no longer a 'heap').
    // You have to either call 'clear()' or 'prepare()' before
    // top/push/pop.

    const _Sequence& sort() {


      sort_heap(c.begin(),c.end(),comp);

      is_valid = false;

      return c;

    }

    void setValueLimit(const limit_type& _new_max_val) {

      if( empty() || limit_comp(top(),_new_max_val) )
	max_val_limit = _new_max_val;
      else {
	dbg::assertion(dbg::error, DBG_ASSERTION( false ));
	throw BoundPriorityQueueError("New value limit is already exceeded by existing elements.");
      }
    }

    void setSizeLimit(size_type _new_size_limit) {

      if( c.size() <= _new_size_limit )
	size_limit = _new_size_limit;
      else {
	dbg::assertion(dbg::error, DBG_ASSERTION( false ));
	throw BoundPriorityQueueError("New size limit is already exceeded by the container.");
      }
      ATLOG_ASSERT_1(size_limit > 0);
    }


    size_type getSizeLimit() const {
      return size_limit;
    }

    void 
    push(const value_type& x) 
    {

      if(! limit_comp(x,max_val_limit))
	return;

      dbg::assertion(dbg::error, DBG_ASSERTION( is_valid ));

      if( ! is_valid )
	throw BoundPriorityQueueError("Invalid state for top/push/pop operations");

      if(c.size() >= size_limit) {
      
	if( size_limit > 0 && comp(x,top()) )
	  pop();
	else
	  return;
      }

      try 
	{
	  c.push_back(x); 
	  push_heap(c.begin(), c.end(), comp);
	}
      catch(...)
	{
	  c.clear();
	  throw; 
	}
    }

    // This function is defined so that back_inserter(...) would work

    void push_back(const value_type& x) {

      push(x);

    }

    void 
    pop() 
    {

      dbg::assertion(dbg::error, DBG_ASSERTION( is_valid ));

      if( ! is_valid )
	throw BoundPriorityQueueError("Invalid state for push/pop operations");
      try 
	{
	  pop_heap(c.begin(), c.end(), comp);
	  c.pop_back();
	}
      catch(...)
	{
	  c.clear();
	  throw; 
	}
    }
  }; // class BoundPriorityQueue


} // namespace PRODDL


#endif // PRODDL_COMMON_QUEUE_H__
