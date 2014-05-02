//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_ASA_ADAPTORS_H__
#define AT_ASA_ADAPTORS_H__

//interface classes to ASA optimization package
//based on asa_adaptor

#include "asa_adaptor.hpp"

#include "bzalgor.h"
#include "blitzdef.h"
#include "contdef.h"
#include "gridcut.h"
#include "eulerang.h"
#include <boost/timer.hpp>

inline void bracket_vectors(const Vvect& vr, vect& vmin, vect& vmax)
{
  using namespace blitz_ext;
  if(vr.size() > 0) {
    vmin = vr(0);
    vmax = vr(0);
    for(int i=1; i<vr.size(); i++)
      bracket(vr(i),vmin,vmax);
  }
}

//find min and maximum distance between different points in an array:
inline void bracket_dist(const Vvect& vr, mfl& dmin, mfl& dmax)
{
  if(vr.size() >= 2) {
    dmin = mods(vect(vr(0) - vr(1)));
    dmax = dmin;
    for(int i=2 ; i < vr.size(); i++)
      for(int j=0; j < i; j++)
	{
	  mfl d2 = mods(vect(vr(i) - vr(j)));
	  if(d2 < dmin) dmin = d2;
	  else if(d2 > dmax) dmax = d2;
	}
    dmin = sqrt(dmin);
    dmax = sqrt(dmax);
  }
}

//"two rigid bodies" optimization - parameters are
//xyz coords of center of mass and three euler angles.
//cost function is atom-atom potential

template<class Fcost> class two_body_asa : public asa_adaptor<two_body_asa<Fcost> >
{
 protected:
  enum { nparam = 6, nbody = 2 };
  //object that calculates energy of interaction between two points:
  Fcost fe;
  //energy calculating object prefilled with coordinates of
  //first body in constructor:
  GridCutoff gridcut;
  //initial position of the center of mass of both bodies:
  vvect vcm;
  //min and max coordinates of both bodies:
  vvect vmin;
  vvect vmax;
  //maximum distance between points in each body:
  fvect fdmax;
  //initial coordinates of both bodies:
  vVvect vr;
  //buffer for second body coordinates
  Vvect vr2;
  //box where make sense to move the center of mass of the second body,
  //taking into account the ultimate goal of locating the position
  //of maximum contact between the bodies
  vect cm2min, cm2max;
  //used to break computation on time limit
  boost::timer timer_opt;
  //time limit in sec:
  double time_opt_max;
  //harmonic potential used to bring bodies together:
  struct cm_harm
  {
    //max value
    mfl max_val;
    mfl max_dist; //max distance
    mfl k; //harmonic force constant
    cm_harm() : max_dist(0), k(0) { }
    cm_harm(mfl _max_dist, mfl _max_val=1.0) :
      max_val(_max_val), max_dist(_max_dist)
      {
	if(max_dist != 0) k = max_val/max_dist;
	else k = 0;
      }
    mfl f(mfl dist) {
      return k * dist;
    }
  };
  cm_harm sqeezer;
  //if calculated cost function value is greater than CostConstrMax,
  //we declare that value invalid in f().
  //Note: must not be less than sqeezer.max_val
  mfl CostConstrMax;
public:
  typedef asa_adaptor<two_body_asa> Base;
  two_body_asa(const Vvect& _v1, const Vvect& _v2, const Fcost& _fe,  mfl cutoff = 9.0) : 
    Base(this), fe(_fe), vcm(nbody), vmin(nbody), vmax(nbody),
    fdmax(nbody),
    vr(nbody), vr2(_v2.size()), time_opt_max(fHUGE)
    {
      asa_param.resize(nparam);
      asa_lower_bound.resize(nparam);
      asa_upper_bound.resize(nparam);
      //DEBUG: use these parameter settings with ASA_TEST cost function
      //in order to reproduce original ASA results. Also, change nparam to 4.
      //asa_lower_bound = -10000;
      //asa_upper_bound = 10000;
      //asa_param = 999,-1007,1001,-903;
      vr[0].reference(_v1.copy());
      vr[1].reference(_v2.copy());
      for(int i = 0; i < nbody; i++) {
	vcm[i] = mean(vr[i]);
	bracket_vectors(vr[i],vmin[i],vmax[i]);
	mfl dummy;
	bracket_dist(vr[i],dummy,fdmax[i]);
      }
      //(vmin[0], vmax[0]) gives us the enclosing box for the first body,
      //which is kept still during the optimization. Then the smallest
      //enclosing box for both bodies during the optimization is:
      vect vboxmin = vmin[0] - fdmax[1], vboxmax = vmax[0] + fdmax[1];
      //initialize grid parameters:
      gridcut.reset(vboxmin,vboxmax,cutoff,vr[0].size());
      //fill grid with coordinates of the first body:
      gridcut.Ef(vr[0],true,fe);
      //looks like the safe guess and must give some room for
      //not-so-close contacts:
      cm2min = vboxmin*5; cm2max = vboxmax*5;
      //TODO: to get the strict estimate of max distance between
      //centers of mass with the established boundary box,
      //we need to find the max distance to the vertices from first body's
      //cm. The simpler way below gives the upper estimate:
      sqeezer = cm_harm(mod(cm2max - cm2min));
      CostConstrMax = sqeezer.max_val*10;
    }
  //basic cost function:
  double f(const vect& angles, const vect& xyz) {
    //put new atom coords into vr2 buffer:
    move_r(angles,xyz,vr2);
    gridcut.reset_n_out();
    mfl f_nb = gridcut.Ef(vr2,false,fe);
    //add contraction potential in order to
    //prevent ASA from stucking at plato:
    vect vcm2 = vcm[1] + xyz; //new cm position
    mfl d02 = mod(vect(vcm[0]-vcm2));
    mfl f_contr = sqeezer.f(d02);
    mfl f_cost = f_nb + f_contr;
    return f_cost;
  }
  double f(double *x,
	   double *parameter_lower_bound,
	   double *parameter_upper_bound,
	   double *cost_tangents,
	   double *cost_curvature,
	   ALLOC_INT * parameter_dimension,
	   int *parameter_int_real,
	   int *cost_flag,
	   int *exit_code,
	   USER_DEFINES * USER_OPTIONS) {
    //convert parameters to rigid body coords:
    vect angles(x[0],x[1],x[2]);
    vect xyz(x[3],x[4],x[5]);
    mfl f_cost = f(angles,xyz);
    //ASA bails out if cost is too large (1e18 by default)
    if(f_cost > MAX_DOUBLE) f_cost = MAX_DOUBLE-10.0;
    if(f_cost < -MAX_DOUBLE) f_cost = -MAX_DOUBLE+10.0;
    //impose constraints on parameter space by setting
	//const function as having invalid value - may be
	//better than climbing high hills
    if(f_cost > CostConstrMax) *cost_flag = FALSE;
    else *cost_flag = TRUE;
    if(timer_opt.elapsed() > time_opt_max) {
      //I don't know another way to stop ASA but
      //to cause BAD_COST_FUNCTION error status:
      //f_cost = MAX_DOUBLE*1.1;
      USER_OPTIONS->Immediate_Exit = TRUE;
    }
    //DEBUG:
    //mfl d01 = mod(vect(vcm[0]-vcm[1]));
    //mfl f_ini = gridcut.ENB(vr[1],false);
    //OUTVAR(f_cost); OUTVAR(angles); OUTVAR(xyz);
    //OUTVAR(vcm[0]); OUTVAR(vcm[1]); OUTVAR(vcm2);
    //OUTVAR(d01); OUTVAR(d02); OUTVAR(f_ini); 
    //OUTVAR(f_nb); OUTVAR(f_contr); cout << endl;
    return f_cost;
  }
  //takes euler angles and cm position (like in gramm, xyz is
  //the displacment of cm from the original structure, 
  //so you can call run(angles,xyz)
  //if xyz is taken from gramm .res file.
  //max_time_minute - limit on computation time in minutes
  void run(vect& angles, vect& xyz, double max_time_minute=fHUGE/1000) {
    const mfl pi = M_PI, pi2 = 2 * M_PI;
    timer_opt.restart();
    time_opt_max = max_time_minute*60;
    asa_param = angles[0], angles[1], angles[2],
                xyz[0], xyz[1], xyz[2];
    vect xyzmin = cm2min - vcm[1];
    vect xyzmax = cm2max - vcm[1];
    asa_lower_bound = -pi2, -pi, -pi2, 
                      xyzmin[0], xyzmin[1], xyzmin[2];
    asa_upper_bound = pi2, pi, pi2,
                      xyzmax[0], xyzmax[1], xyzmax[2];

    run_asa();

    angles = vect(asa_param(0),asa_param(1),asa_param(2));
    xyz    = vect(asa_param(3),asa_param(4),asa_param(5));
  }
  //put into vr_new cartesian coords of the second body corresponding to
  //(angles, xyz) coords
  void move_r(const vect& angles,const vect& xyz, Vvect& vr_new) const {
    matr3 mrot;
    AtGeom::eu_xyz::EuRotMatr(angles,mrot);
    //put new 3D coords into vr:
    ATALWAYS(vr[1].size() == vr_new.size(), "two_body_asa::move_r() - array size mismatch");
    AtGeom::EulerMove(mrot,xyz,vcm[1],vr[1],vr_new);    
  }
};

#endif //AT_ASA_ADAPTORS_H__
