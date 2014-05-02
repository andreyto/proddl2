#include <math.h>

namespace DockTK { namespace nr {

template<typename T_num>
T_num betacf(T_num a, T_num b, T_num x)
{
  const T_num  EPS = 3.0e-7;
  const T_num  FPMIN = 1.0e-30;
  const int    MAXIT = 100;
	void nrerror(char error_text[]);
	int m,m2;
	T_num aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}

} } // namespace DockTK::nr


// instantiate

namespace DockTK { namespace nr {

template float betacf(float , float , float );
template double betacf(double , double , double );

} } // namespace DockTK::nr

/* (C) Copr. 1986-92 Numerical Recipes Software =$jX16#S'j#(+1. */
