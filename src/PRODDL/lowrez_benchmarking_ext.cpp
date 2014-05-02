//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include <iterator>

#include "PRODDL/External/nr/nr.hpp"

#include <math.h>

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/to_string.hpp"

#include "PRODDL/Common/math.hpp"

#include "PRODDL/Common/bz_vect_ext.hpp"

#include "PRODDL/Common/common_types.hpp"

#include "PRODDL/Pyarray/array.hpp"

#include "PRODDL/Blitz/bzarray_iter.hpp"

#include "PRODDL/Common/iterator.hpp"

#include "PRODDL/Common/function.hpp"


namespace PRODDL {

template<typename T_num> struct point_traits {
  typedef typename PRODDL::common_types::point_type<T_num,3>::Type Type;
};

template<typename T_num>
T_num BinomialCPF(int n,int k,T_num p)
{
  //special cases:
  ATALWAYS( !(n<0 || k<0 || p<0),std::string("Negative arguments: k=")+PRODDL::to_string(k)+" p="+PRODDL::to_string(p)+" n="+PRODDL::to_string(n));
  if(k==0) return 1;
  if(n==0) return 0;
  if(p==0) return 0;
  //general case:
  return PRODDL::nr::betai(T_num(k),T_num(n-k+1),p);  
}


template<typename T_num>
T_num RndHitSphere(T_num R,T_num rsite,int n,int k)
{
  //check range of arguments for validity:
  ATALWAYS(R>0,"R="+PRODDL::to_string(R));
  ATALWAYS(rsite>=0, "rsite="+PRODDL::to_string(rsite));
  ATALWAYS(n>=0,"n="+PRODDL::to_string(n));
  ATALWAYS(k>=0,"k="+PRODDL::to_string(k));
  //the line below is exact value of
  //the cap surface on the sphere.
  T_num site_surf = M_PI*PRODDL::Math::pow2(rsite);
  T_num tot_surf  = 4*M_PI*PRODDL::Math::pow2(R);
  T_num p = site_surf/tot_surf;
  return BinomialCPF(n,k,p);
}


//Count weight of hits of xyz vectors inside sphere of radius rS
//around site at position vS. [f_vr,l_vr) - xyz positions
//[f_w,f_w+(l_vr-f_vr)) - statistical weights, rS - radius of Site,
//vS - position of Site. 
template<class InpIter1,class InpIter2,typename T_num> 
typename std::iterator_traits<InpIter2>::value_type 
CntHitSph(InpIter1 f_vr,InpIter1 l_vr,InpIter2 f_w,T_num rS,
const typename point_traits<T_num>::Type& vS)
{
  T_num rS2 = PRODDL::Math::pow2(rS);
  typename std::iterator_traits<InpIter2>::value_type whit = 0;
  for( ; f_vr != l_vr; ++f_vr, ++f_w)
    {
      if(blitz_ext::dotSelf(*f_vr-vS) <= rS2) whit += *f_w;
    }
  return whit;
}


//check for 'good docking' by original criterion of
//'random hits on a sphere' model. [f_vr,l_vr) - xyz positions
//[f_w,f_w+(l_vr-f_vr)) - statistical weights, rS - radius of Site,
//vS - position of S., sum_w - sum(f_w...), transfered here from outside
//so that this function can be more efficiently used inside some
//inner loop where just vS changes; p_val - 'pi value'; rRND - radius
//of sphere which carries random matches.
//Results: w_hit -
//number (weight) of hits; p_rand - probability to get n_hit by chance;
//Return: true if n_hit is statisticaly significant.
//InpIterW::value_type must be int, since Bernully distribution used
//in RndHitSphere does not make sence otherwise.
template<class InpIterV,class InpIterW, typename T_num> 
bool SignifRndSph(InpIterV f_vr,InpIterV l_vr,InpIterW f_w,
		  T_num rS,
		  const typename point_traits<T_num>::Type& vS,
		  T_num p_val,
		  T_num rRND,
		  int sum_w,int& w_hit,
		  T_num& p_rand)
{
  w_hit = CntHitSph(f_vr,l_vr,f_w,rS,vS);
  p_rand = RndHitSphere(rRND,rS,sum_w,w_hit);
  return p_rand < p_val;
}

} // namespace PRODDL {


