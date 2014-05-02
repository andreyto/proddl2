//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_ITERATOR_H__
# define AT_ITERATOR_H__

#include <iterator>

namespace PRODDL {

  //"Pipe" iterator: generalization
  //of "ostream_iterator" usage:
  //Suppose that we have some unary function object of type F, 
  //i.e. the one that has operator() (const Y& y)
  //working on some other object of type Y.
  //Suppose, also, that some algorithm A produces
  //objects of type Y and copies them into an output range
  //represented by output iterator OutputIterator (e.g.
  //"std::transform(...)"). Then in places where we need to
  //write "A(...,OutputIterator result)", we can write
  //"A(...,piper(f))", where f is object of type F and
  //"piper()" is a helper function creating pipe_iterator<F>(f).
  //In other words, if p is pipe_iterator<F>(f), then f(y) is equal to
  //*p++=y.
  //Why we need it: same as with "pipe" in general:
  //if A is parametrized by OutputIterator type, 
  //then we have a freedom of either storing output of
  //A into some container, or directly feeding A's output to
  //F without intermediate storage.
  template <class UnaryFunction>
  class opipe_iterator 
  {
  protected:
    UnaryFunction* function;
  public:
    typedef UnaryFunction       function_type;
    typedef std::output_iterator_tag iterator_category;
    typedef void                value_type;
    typedef void                difference_type;
    typedef void                pointer;
    typedef void                reference;

    explicit opipe_iterator(UnaryFunction& x_) : function(&x_) {}
    template<class T> opipe_iterator<UnaryFunction>&
    operator=(const T& value_) 
    { 
      (*function)(value_);
      return *this;
    }
    opipe_iterator<UnaryFunction>& operator*() { return *this; }
    opipe_iterator<UnaryFunction>& operator++() { return *this; }
    opipe_iterator<UnaryFunction>& operator++(int) { return *this; }
  };

  template <class UnaryFunction>
  inline opipe_iterator<UnaryFunction> make_opipe_iterator(UnaryFunction& x) 
  {
    return opipe_iterator<UnaryFunction>(x);
  }


  //input pipe iterator: if g is generator object,
  //i is ipipe_iterator<generator>(g,...) and
  //y is of some other type, then y = *i++ is equal
  //to: y = g(). Size_t argument to constructor is
  //used to create a counter, incremented with each
  //++ operation, and test two iterators for equality.
  //Two iterators are equal if their counters are equal.
  //Essentialy, it's just a counter inside operator++(),
  //combined with generator invocation inside operator*(),
  //and rules for InputIterator concept must be strictly
  //followed when using this type.

  template <class _Generator>
  class ipipe_iterator 
  {
  public:
    typedef _Generator          function_type;
    typedef std::input_iterator_tag  iterator_category;
    typedef typename function_type::result_type value_type;
    typedef int                 difference_type;
    typedef const value_type*   pointer;
    typedef const value_type&   reference;

  protected:
    _Generator* function;
    size_t i;
  public:
    explicit ipipe_iterator(_Generator& _x, size_t _i=0) : 
      function(&_x), i(_i) {}
    ipipe_iterator(const ipipe_iterator& x) :
      function(x.function), i(x.i)
    {
    }
    ipipe_iterator& operator= (const ipipe_iterator& x)
    {
      function = x.function;
      i = x.i;
      return *this;
    }
    value_type operator*() 
    { 
      return (*function)(); 
    }
    pointer operator->() { return &(operator*()); }  
    ipipe_iterator& operator++() { 
      i++;
      return *this;
    }
    ipipe_iterator operator++(int)  {
      ipipe_iterator __tmp(*this);
      i++;
      return __tmp;
    }
    bool equal(const ipipe_iterator& x) const
    {
      return (i == x.i);
    }
  };

  template <class _Generator> inline bool operator==
  (const ipipe_iterator<_Generator>& x1,const ipipe_iterator<_Generator>& x2)
  {
    return x1.equal(x2);
  }

  template <class _Generator> inline bool operator!=
  (const ipipe_iterator<_Generator>& x1,const ipipe_iterator<_Generator>& x2)
  {
    return !x1.equal(x2);
  }


  template <class _Generator>
  inline ipipe_iterator<_Generator> make_ipipe_iterator(_Generator& x,size_t index=0) 
  {
    return ipipe_iterator<_Generator>(x,index);
  }

} // namespace PRODDL

#endif // AT_ITERATOR_H__
