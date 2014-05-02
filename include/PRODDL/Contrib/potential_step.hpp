//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

  class MultiStep
  {

  public:

    typedef typename Geom::SpaceTraits<T_num>::Point3 Point;

  protected:

    // raw parameters

    T_num sigma; // standard LJ well radius parameter
    T_num eps; // standard LJ barrier height parameter
    T_num alpha; // smoothing parameter
    T_num cut; // cutoff
    T_num cutP2; // cutoff^2

    // derived values constant between calls to value and gradient methods

    T_num _s6;
    T_num _as6;
    T_num _4es6;
    T_num _4es12;
    T_num _12es6;

    // shifted energy

    T_num _fshift;

  public:
    
    SoftCoreLJ() 
    {}

    SoftCoreLJ(T_num _sigma, T_num _eps, T_num _alpha, T_num _cut):
      sigma(_sigma),
      eps(_eps),
      alpha(_alpha),
      cut(_cut),
      cutP2(_cut*_cut)
    {
      _s6 = Math::pow6(sigma);
      _as6 = alpha*_s6;
      _4es6 = 4.*eps*_s6;
      _4es12 = _4es6*_s6;
      _12es6 = 3.*_4es6;

      _fshift = 0;
      _fshift = f2(cutP2);
    }

    // Get square of the distance and return the value of the potential
    // Requirements: assumes that rP2 is always less than xP2 (see projection.hpp
    // for projection function that will call this f2() method after checking
    // that rP2 is less than value returned by this object's getRadius() method.

    T_num f2(T_num rP2) const {
      T_num A = Math::pow3(rP2) + _as6;
      return (_4es12/A - _4es6)/A - _fshift;
      // DEBUG:
      //return rP2;
    }

    Point g2(T_num rP2,const Point& xyz) const {
      T_num _r6 = Math::pow3(rP2);
      T_num A = _r6 + _as6;
      T_num B = _s6/(_r6 + _as6) - 1.;
      T_num C = -6.*Math::pow2(rP2)*(_4es6*A*B + _4es12)/Math::pow3(A);
      // DEBUG:
      return xyz*C;
      //return 2*xyz;
    }

    // Return the radius beyond which this potential is always zero.

    T_num getRadius() const {
      return cut;
    }

    // return getRadius()^2
    
    T_num getRadiusP2() const {
      return cutP2;
    }

  };





T_num customEnergy(T_num alpha, T_num sigma, T_num epsilon, T_num delta, T_num r12)
{
  /* alpha, sigma, epsilon are the parameters of a LJ potential */
  /* delta is a step of a stair-wise potential*/

  

  T_num l,r1,r2,v1,v2,v3, energy;

  l=std::floor(r12/delta);

  r1=l*delta;
  r2=(l+1.0)*delta;

  v1=Math::pow3(r2)/(Math::pow6(r2)+alpha*Math::pow6(sigma))-Math::pow3(r1)/(Math::pow6(r1)+alpha*Math::pow6(sigma));

  v2=std::atan(Math::pow3(r2)/(std::sqrt(alpha)*Math::pow3(sigma)))-std::atan(Math::pow3(r1)/(std::sqrt(alpha)*Math::pow3(sigma)));

  v3=4.0*epsilon*Math::pow3(sigma)/(Math::pow3(r2)-Math::pow3(r1));

  energy=v3*(v1*Math::pow3(sigma)/(2.0*alpha)+v2*(-1.0+1.0/(2.0*alpha))/std::sqrt(alpha));

  return energy;
}

