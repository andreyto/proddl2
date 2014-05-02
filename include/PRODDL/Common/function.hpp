//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_FUNCTION_H__
# define AT_FUNCTION_H__

#include <stdexcept>

namespace PRODDL {

  // Generator classes

  // base class for generators

  template<class _ResultType> struct generator
  {
    typedef _ResultType result_type;
  };


  //generator that generates same constant value

  template<class T> struct gen_const : public generator<T>
  {
    T x;
    gen_const(const T& _x) : x(_x) {}
    const T& operator() () { return x; }
  };

  //helper function to obtain genconst object

  template<class T> gen_const<T> make_gen_const(const T& _x)
  {
    return gen_const<T>(_x);
  }


  // 
//   template<class T_map> 
//   typename T_map::data_type& 
//   getMapValue(T_map& m, const typename T_map::key_type& key) throw(std::out_of_range)
//   {
//     typename T_map::iterator iter = m.find(key);
//     if( iter == m.end() )
//       throw std::out_of_range("getMapValue(): key not found in map.");
//     return (*iter).second;
//   }

//   template<class T_map> 
//   const typename T_map::data_type& 
//   getMapValue(const T_map& m, const typename T_map::key_type& key)  throw(std::out_of_range)
//   {
//     return getMapValue(const_cast<T_map&>(m),key);
//   }  

  template<class T_map, class KeyType> 
  typename T_map::data_type& 
  getMapValue(T_map& m, const KeyType& key) throw(std::out_of_range)
  {
    typename T_map::iterator iter = m.find(key);
    if( iter == m.end() )
      throw std::out_of_range("getMapValue(): key not found in map.");
    return (*iter).second;
  }

  template<class T_map, class KeyType> 
  const typename T_map::data_type& 
  getMapValue(const T_map& m, const KeyType& key)  throw(std::out_of_range)
  {
    return getMapValue(const_cast<T_map&>(m),key);
  }  


} // namespace PRODDL

#endif // AT_FUNCTION_H__
