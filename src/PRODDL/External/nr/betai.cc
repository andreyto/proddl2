#include <math.h>

#include "DockTK/External/nr/nr.hpp"

namespace DockTK { namespace nr {

template<typename T_num>
T_num betai(T_num a, T_num b, T_num x)
{
  //T_num betacf(T_num a, T_num b, T_num x);
  //T_num gammln(T_num xx);
	void nrerror(char error_text[]);
	T_num bt;

	if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,T_num(1.0)-x)/b;
}

} } // namespace DockTK::nr


// instantiate

namespace DockTK { namespace nr {

template float betai(float , float , float );
template double betai(double , double , double );

} } // namespace DockTK::nr

/* (C) Copr. 1986-92 Numerical Recipes Software =$jX16#S'j#(+1. */
