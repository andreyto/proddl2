#include <math.h>

#include <limits>

#include "DockTK/External/nr/nrutil.hpp"

namespace DockTK { namespace nr {


  template<typename T_num> struct ZbrentConstants
  {
    static T_num EPS() throw() {  return std::numeric_limits<T_num>::epsilon(); }
  };

  template<> struct ZbrentConstants<float>
  {
    static float EPS() throw() { return 3.0e-8; }
  };


template<typename T_num>
T_num zbrent(T_num (*func)(T_num), T_num x1, T_num x2, T_num tol)
{
  const int ITMAX = 100;
  const T_num  EPS = ZbrentConstants<T_num>::EPS();
	int iter;
	T_num a=x1,b=x2,c=x2,d=0,e=0,min1,min2;
	T_num fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		nrerror("Root must be bracketed in zbrent");
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(b);
	}
	nrerror("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}

} } // namespace DockTK::nr


// instantiate

namespace DockTK { namespace nr {

  template float zbrent(float (*func)(float), float x1, float x2, float tol);
  template double zbrent(double (*func)(double), double x1, double x2, double tol);

} } // namespace DockTK::nr


/* (C) Copr. 1986-92 Numerical Recipes Software =$jX16#S'j#(+1. */
