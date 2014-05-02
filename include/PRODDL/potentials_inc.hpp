//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PRODDL_POTENTIALS_H__
#  error POTENTIALS_INC.HPP MUST BE INCLUDED FROM WITHIN POTENTIALS.HPP
#endif

#ifndef AT_PRODDL_POTENTIALS_INC_H__
#define AT_PRODDL_POTENTIALS_INC_H__

// "Radial potential" concept.
//
// Radial potentials must conform to the following model:
// Potential is defined by the distance to the point in space where
// we want to compute the potential.
// Each potential has a radius beyond which the potential is zero (cutoff radius).
//
// Methods:
//
// T_num getRadius() must return the radius of the potential.
//
// T_num getRadiusP2() must return getRadius()^2 (we assume that this object will cache this value inside itself)
//
// T_num f2(T_num rP2) takess the square of the distance 'rP2' and returns the value of the potential
// Requirements: then this method is called it must be always guaranteed that rP2 is less than getRadius()^2.
// It is the responsibility of the caller to compare rP2 with the square of the potential's radius and
// call 'f2()' only if rP2 <= getRadius()^2.
// Rational: This method is designed to be called from various loops accumulating pairwise potential interaction between
// a large number of points (see "grid/project.hpp" for the examples of functions which accept such objects).
//


// "Radial potential with 1st derivative" concept.
// Conforms to the "Radial potential" concept plus provides a 1st derivative (gradient)
//
// Methods:
//
// Point g2(T_num rP2,const Point& xyz) takes the square of a distance 'rP2' and the vector xyz, and returns
// the gradient.
// Requirements: same as for 'f2()' method, plus 'xyz^2 == rP2'



// 'Step function' : 'y' inside 'x', zero - elsewhere
// Concept: "Radial potential"

class Step
{

protected:
	T_num x;
	T_num xP2;
	T_num y;

public:

	Step(T_num _x, T_num _y):
	  x(_y),
		  xP2(_x*_x),
		  y(_y)
	  {
	  }

	  // Get square of the distance and return the value of the potential
	  // Requirements: assumes that rP2 is always less than xP2 (see projection.hpp
	  // for projection function that will call this f2() method after checking
	  // that rP2 is less than value returned by this object's getRadius() method.

	  T_num f2(T_num rP2) {
		  return y;
	  }

	  // Return the radius beyond which this potential is always zero.

	  T_num getRadius() const {
		  return x;
	  }

	  // return getRadius()^2

	  T_num getRadiusP2() const {
		  return xP2;
	  }

};


// 'TwoSteps function' : 'y1' inside 'x1', 'y2' between 'x1' and 'x2', zero - elsewhere
// Concept: "Radial potential"

class TwoSteps
{

protected:
	Step step1, step2;

public:

	TwoSteps(T_num _x1, T_num _y1, T_num _x2, T_num _y2):
	  step1(_x1,_y1), step2(_x2,_y2)
	  {
	  }

	  // Get square of the distance and return the value of the potential
	  // Requirements: assumes that rP2 is always less than x2^2

	  T_num f2(T_num rP2) {
		  if( rP2 <= step1.getRadius() )
			  return step1.f2(rP2);
		  else 
			  return step2.f2(rP2);
	  }

	  // Return the radius beyond which this potential is always zero.

	  T_num getRadius() const {
		  return step2.getRadius();
	  }

};


enum { CUT_POLICY_NONE, CUT_POLICY_SHIFT_ENERGY, CUT_POLICY_SHIFT_FORCE };

// How LJ sigma values for a pair of MM force field atom types are mixed,
// according to either Amber or OPLS default rules, as in 
// the third value returned by MMTK.ForceFields.AmberData.ljParameters(),
// and used by MMTK.ForceFields.NonBondedInteractions.LJForceField.
// 0 means arithmetic average (Amber), and 1 means geometric average (OPLS)

enum { LJ_MIX_0 = 0, LJ_MIX_1 = 1 };

// 'Softcore Lennard Jones' : 6-12 LJ with additional smoothing parameter
// to make a softer repulsion part and exclude discontinuity at 0
// Concept: "Radial potential with 1st derivative"

class SoftCoreLJ
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

	// shifted force

	T_num _gshift;

public:

	SoftCoreLJ() 
	{}

	SoftCoreLJ(T_num _sigma, T_num _eps, T_num _alpha, T_num _cut, int cutPolicy=CUT_POLICY_SHIFT_FORCE):
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
		// Not a mistake: _fshift == 0 is used inside a call to f2() in the next line
		_fshift = 0.;
		_gshift = 0.;
		if(cutPolicy != CUT_POLICY_NONE) {
			Point _gcut;
			_gcut = 0;
			_gcut(0) = cut;
			_gcut = g2(cutP2,_gcut);
			_fshift = - f2(cutP2);
			if(cutPolicy == CUT_POLICY_SHIFT_FORCE) {
				_gshift = - _gcut(0);
				_fshift -= _gshift * cut;
			}
		}
	}

	// Get square of the distance and return the value of the potential
	// Requirements: assumes that rP2 is always less than xP2 (see projection.hpp
	// for projection function that will call this f2() method after checking
	// that rP2 is less than value returned by this object's getRadius() method.

	T_num f2(T_num rP2) const {
		T_num A = Math::pow3(rP2) + _as6;
		return (_4es12/A - _4es6)/A + _fshift + _gshift * std::sqrt(rP2);
		// DEBUG:
		//return rP2;
	}

	Point g2(T_num rP2,const Point& xyz) const {
		T_num _r6 = Math::pow3(rP2);
		T_num A = _r6 + _as6;
		T_num B = _s6/(_r6 + _as6) - 1.;
		T_num C = -6.*Math::pow2(rP2)*(_4es6*A*B + _4es12)/Math::pow3(A);
		// DEBUG:
		return xyz*C + xyz*(_gshift/std::sqrt(rP2)) ;
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


// Repulsive part of soft core LJ potential - for FFT projection

class SoftCoreLJRep
{

public:

	typedef typename Geom::SpaceTraits<T_num>::Point3 Point;

protected:

	// raw parameters

	T_num sigma; // standard LJ well radius parameter
	T_num sigmaAlpha; // averaged sigma to be combined with alpha
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


public:

	SoftCoreLJRep() 
	{}


	// Difference with regular SoftCoreLJ: _sigma must be sqrt(sigma_i), same for _eps,
	// but _sigmaAlpha must be full average sigma.
	// Values of this potential projected onto receptor grid will be multiplied
	// by the values from the ligand grid, which will be: sqrt(sigma_j)^12 * sqrt(eps_j).
	// The combined potential RepR * RepL + AttrR * AttrL will give an exact value
	// of a pairwise LJ potential assuming OPLS default rule of combining sigma values
	// for pairs of atom types (sqrt(sigma_i*sigma_j)) and alpha == 0.
	// When alpha > 0, the difference from the exact pairwise potential will still be small.

	SoftCoreLJRep(T_num _sigma, T_num _sigmaAlpha, T_num _eps, T_num _alpha, T_num _cut):
	sigma(_sigma),
		sigmaAlpha(_sigmaAlpha),
		eps(_eps),
		alpha(_alpha),
		cut(_cut),
		cutP2(_cut*_cut)
	{
		_s6 = Math::pow6(sigma);
		_as6 = alpha*Math::pow6(sigmaAlpha);
		_4es6 = 4.*eps*_s6;
		_4es12 = _4es6*_s6;
		_12es6 = 3.*_4es6;


	}

	// Get square of the distance and return the value of the potential
	// Requirements: assumes that rP2 is always less than xP2 (see projection.hpp
	// for projection function that will call this f2() method after checking
	// that rP2 is less than value returned by this object's getRadius() method.

	T_num f2(T_num rP2) const {
		T_num A = Math::pow3(rP2) + _as6;
		return _4es12/A/A;
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


// Attractive part of soft core LJ, to be used in FFT projection

class SoftCoreLJAttr : public SoftCoreLJRep
{

public:


	SoftCoreLJAttr() 
	{}

	// See notes for the parent class ctor

	SoftCoreLJAttr(T_num _sigma, T_num _sigmaAlpha, T_num _eps, T_num _alpha, T_num _cut) :
	SoftCoreLJRep(_sigma,_sigmaAlpha,_eps,_alpha,_cut)
	{}


	T_num f2(T_num rP2) const {
		T_num A = Math::pow3(rP2) + SoftCoreLJRep::_as6;
		return -SoftCoreLJRep::_4es6/A;
	}  

};


struct PotAce {

	T_num b, k, a;

	T_num cut;

	// shifted energy

	T_num _fshift;

	// shifted force

	T_num _gshift;


	PotAce() {}

	PotAce(T_num _b, T_num _k, T_num _a, T_num _cutoff, int cutPolicy=CUT_POLICY_SHIFT_FORCE) :
	b(_b), k(_k), a(_a), cut(_cutoff)
	{
		// Not a mistake: _fshift == 0 is used inside a call to f2() in the next line
		_fshift = 0.;
		_gshift = 0.;
		if(cutPolicy != CUT_POLICY_NONE) {
			T_num cutP2 = cut*cut;
			Point _gcut;
			_gcut = 0;
			_gcut(0) = cut;
			_gcut = g2(cutP2,_gcut);
			_fshift = - f2(cutP2);
			if(cutPolicy == CUT_POLICY_SHIFT_FORCE) {
				_gshift = - _gcut(0);
				_fshift -= _gshift * cut;
			}
		}
	}

	T_num f2(T_num rP2) {
		T_num r = std::sqrt(rP2);
		return 1./(1 + std::exp(b*(k*r-a))) + _fshift + _gshift * r;
	}

	Point g2(T_num rP2, const Point& xyz) {
		T_num r = std::sqrt(rP2);
		T_num bkr = std::exp(b*(k*r-a));
		T_num d = - b*k*bkr/Math::pow2(1+bkr);
		return xyz*(d/r) + xyz*(_gshift/r);
	}

};


// Harmonic potential - used to create a shallow basin
// in order to confine a system to a particular area in space
// Concept: "Infinite radial potential with 1st derivative"

class Harmonic
{

public:

	typedef typename Geom::SpaceTraits<T_num>::Point3 Point;

protected:

	// raw parameters

	T_num k; // harmonic force constant

public:

	Harmonic() 
	{}

	Harmonic(T_num _k):
	k(_k)
	{
	}

	// Get square of the distance and return the value of the potential
	T_num f2(T_num rP2) {
		return k*rP2;
	}

	Point g2(const Point& xyz) {
		return 2*k*xyz;
	}

};

#endif // AT_PRODDL_POTENTIALS_INC_H__
