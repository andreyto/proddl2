#include <math.h>

namespace DockTK { namespace nr {

template<typename T_num>
T_num gammln(T_num xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

} } // namespace DockTK::nr


// instantiate

namespace DockTK { namespace nr {

  template float gammln(float xx);
  template double gammln(double xx);

} } // namespace DockTK::nr


/* (C) Copr. 1986-92 Numerical Recipes Software =$jX16#S'j#(+1. */