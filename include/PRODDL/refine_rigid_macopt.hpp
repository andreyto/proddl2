//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PRODDL_REFINE_RIGID_MACOPT_H__
#define AT_PRODDL_REFINE_RIGID_MACOPT_H__

// Classes to refines docking predictions by local or semi-local minimization

#include "PRODDL/potentials.hpp"

#include "PRODDL/Geom/bounding.hpp"

#include "PRODDL/Optim/optim_mac.hpp"

#include "PRODDL/Optim/objective_func_bz.hpp"

#include "PRODDL/Geom/transformation.hpp"

#include "PRODDL/Geom/move.hpp"

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/g_options.hpp"

namespace PRODDL {

  template<typename T_num>
  class Refinement {

  public:

    typedef typename Potentials<T_num>::PotTotalNonBonded PotTotalNB;
    typedef typename Potentials<T_num>::MolForceParams MolForceParams;
    typedef typename Geom::SpaceTraits<T_num>::Point3 Point;
    typedef typename Geom::SpaceTraits<T_num>::VPoint3 Points;
    typedef Geom::Bounding::Box<T_num> BoundingBox;
    typedef Geom::Bounding::Diameter<T_num> BoundingDiameter;
    typedef typename BoundingBox::PointPair PointPair;

    typedef typename common_types::num_vector_type<T_num>::Type fvect;

    typedef typename Geom::TransformationTraits<T_num>::VRotationTranslation VRotationTranslation;
    typedef Geom::Rotation<T_num> Rotation;
    typedef Geom::Translation<T_num> Translation;
    typedef Geom::RotationTranslation<T_num> RotationTranslation;

    // TODO: Move adaptors to optimization routine into separate class
    typedef Optim::ObjectiveFuncContigBzArr<double,Point,Refinement<T_num>,void> OptimAdaptor;


  protected:

    enum { n_coords = 6 };

    PotTotalNB potTotalNB;

    OptimAdaptor optimAdaptor;

    Optim::OptimConjGradMac optimizer;

    // buffers for coordinate conversions

    Points xyzCoords, xyzGrads;

    // (starting coords - center of starting coords)
    
    Points xyzCoordsZero;

    // center of starting coords

    PointPair rbCoordsZero;

  public:

    Refinement() {}

    // NOTE: Array params should be accepted by values, so that we can safely pass
    // here temporary objects returned by view_as_blitz() Python view creation functions
    // and such !!!

    Refinement(const Points recPoints, const Points ligPoints, const MolForceParams& mfParams):
      optimAdaptor(*this,&Refinement::grad),
      optimizer(optimAdaptor,n_coords) { 

      int runTimeLogLevel;

      gOptions.getdefault("logLevel",runTimeLogLevel,ATLOG_LEVEL_1);

      Logger::setRunTimeLevel(runTimeLogLevel);

      ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() \
		     << ATLOGVAR(ATLOG_LEVEL) \
		     << ATLOGVAR(Logger::getRunTimeLevel()) <<"\n");

      const Options& rigOptions = gOptions.getBlock("rigid");

      T_num tolerance;
      rigOptions.getdefault("tolerance",tolerance,1e-5);

      optimizer.setTolerance(tolerance);

      int iterMax;
      rigOptions.getdefault("iterMax",iterMax,200);

      optimizer.setIterMax(iterMax);


      T_num receptorMovePadding;
      rigOptions.getdefault("receptorMovePadding",receptorMovePadding,1.0);

      // select the size of the partiotioning grid as bounding receptor box
      // along the current coordinate axes ('true' as a 2nd argument to
      // BoundingBox ctor)
      // extended in all directions by the diameter of a ligand

      BoundingBox boundingBox(recPoints,true);
      PointPair bounds = boundingBox.getDiagonal();
      BoundingDiameter boundingDiameter(ligPoints);
      T_num ligSize = boundingDiameter.getSize();
      ligSize += receptorMovePadding;
      bounds[0] -= ligSize;
      bounds[1] += ligSize;

      potTotalNB.init(recPoints,mfParams,bounds,rigOptions);

      xyzCoords.reference(Points(ligPoints.size()));
      xyzGrads.reference(Points(xyzCoords.size()));
      xyzCoordsZero.reference(Points(xyzCoords.size()));
    }

    void grad(const Points& rbArrCoords, Points& rbArrGrads) {

      PointPair rbCoords, rbGrads;
      rbCoords(0) = rbArrCoords(0);
      rbCoords(1) = rbArrCoords(1);
      //xyzCoords = xyzCoordsZero;
      //RotationTranslation tr(rbCoords);
      //tr(xyzCoords);
      Geom::Transforms::RigidBody<T_num>::rigidxyz(xyzCoordsZero,rbCoords,xyzCoords);
      potTotalNB.g(xyzCoords,xyzGrads);
      Geom::Transforms::RigidBody<T_num>::xyzGradToRigid(xyzCoords,
							 xyzGrads,
							 rbCoords,
							 rbGrads);
#if ATDEBUG_LEVEL > 8
      {
	T_num e = potTotalNB.f(xyzCoords);
	T_num epsilon = 1e-10;
	T_num normGrads = std::sqrt(blitz::sum(blitz::sum(xyzGrads*xyzGrads)));
	xyzGrads /= normGrads;
	xyzCoords += epsilon*xyzGrads;
	T_num e1 = potTotalNB.f(xyzCoords);
	T_num delta_e = (e1 - e)/epsilon;
	ATOUTVAR(delta_e); ATOUTVAR(normGrads);
	ATOUTVAR(rbCoords); ATOUTVAR(rbGrads); ATOUTVAR(e); ATOUTENDL();
      }
#endif

      rbArrGrads(0) = rbGrads(0);
      rbArrGrads(1) = rbGrads(1);
    }

    T_num f(const Points& rbArrCoords) {

      PointPair rbCoords;
      rbCoords(0) = rbArrCoords(0);
      rbCoords(1) = rbArrCoords(1);
      //xyzCoords = xyzCoordsZero;
      //RotationTranslation tr(rbCoords);
      //tr(xyzCoords);
      Geom::Transforms::RigidBody<T_num>::rigidxyz(xyzCoordsZero,rbCoords,xyzCoords);
      return potTotalNB.f(xyzCoords);
    }





    T_num 
    refineOne(const Points xyzLigStart,RotationTranslation& resultTransform) {
      //ATOUTVAR(xyzLigStart.size()); ATOUTVAR(xyzCoords.size()); ATOUTVAR(xyzCoordsZero.size()); ATOUTENDL();
      fvect mass(xyzLigStart.size());
      mass = 1.0;
      Geom::Transforms::RigidBody<T_num>::standardOrientation(xyzLigStart,mass,xyzCoordsZero,rbCoordsZero);
      PointPair rbCoordsResult(2);
      rbCoordsResult(0) = rbCoordsZero(0);
      rbCoordsResult(1) = rbCoordsZero(1);
#if ATDEBUG_LEVEL > 8
      {
	Geom::Transforms::RigidBody<T_num>::rigidxyz(xyzCoordsZero,rbCoordsZero,xyzCoords);
	xyzCoords -= xyzLigStart;
	T_num badDiff = std::sqrt(blitz::sum(blitz::sum(xyzCoords*xyzCoords))/xyzCoords.size());
	ATOUTVAR(badDiff); ATOUTENDL();
	RotationTranslation backTrans = Geom::Transforms::RigidBody<T_num>::transformation(rbCoordsZero);
	xyzCoords = xyzCoordsZero;
	backTrans(xyzCoords);
	xyzCoords -= xyzLigStart;
	T_num badDiffTr = std::sqrt(blitz::sum(blitz::sum(xyzCoords*xyzCoords))/xyzCoords.size());
	ATOUTVAR(badDiffTr); ATOUTENDL();	
	ATOUTVAR(rbCoordsResult); ATOUTENDL();
	Points rbArrCoords(2), rbArrGrads(2);
	rbArrCoords(0) = rbCoordsResult(0);
	rbArrCoords(1) = rbCoordsResult(1);
	T_num epsilon = 1e-10;
	T_num e0 = f(rbArrCoords);
	T_num e0_1 = potTotalNB.f(xyzLigStart);
	ATOUTVAR(e0); ATOUTVAR(e0_1); ATOUTENDL();
	grad(rbArrCoords,rbArrGrads);
	T_num normGradsRbc = std::sqrt(blitz::sum(blitz::sum(rbArrGrads*rbArrGrads)));
	Points unitGrads(2);
	rbArrGrads /= normGradsRbc;
	rbArrCoords += epsilon*rbArrGrads;
	T_num e1 = f(rbArrCoords);
	T_num delta_e_rbc = (e1 - e0)/epsilon;
	ATOUTVAR(rbArrGrads); 
	ATOUTVAR(delta_e_rbc); ATOUTVAR(normGradsRbc); ATOUTENDL();
      }
#endif
      optimizer.optimize(reinterpret_cast<T_num*>(&rbCoordsResult),n_coords);
      // make the result to be relative to the initial position
      resultTransform = Geom::Transforms::RigidBody<T_num>::transformation(rbCoordsResult)*
	Geom::Transforms::RigidBody<T_num>::transformation(rbCoordsZero).inverse();
      xyzCoords = xyzLigStart;
      resultTransform(xyzCoords);
      T_num e = potTotalNB.f(xyzCoords);
#if ATDEBUG_LEVEL > 8
      {
	Points rbArrCoords(2);
	rbArrCoords(0) = rbCoordsResult(0);
	rbArrCoords(1) = rbCoordsResult(1);
	T_num e_rb = f(rbArrCoords);
	ATOUTVAR(rbCoordsResult); ATOUTVAR(e); ATOUTVAR(e_rb); ATOUTENDL();
      }
#endif
      return e;
    }

    
    void
    setPositionsReceptor(const Points points) {

      potTotalNB.setPoints1(points);

    }

    void
    refineMany(const Points xyzLigReference, 
	       const VRotationTranslation transforms,
	       VRotationTranslation transforms_new,
	       fvect e_new) {

      Points xyzLigStart(xyzLigReference.size());

      for(int i_trans = 0; i_trans < transforms.size(); i++) {

	xyzLigStart = xyzLigReference;
	const RotationTranslation& transform = transforms(i_trans);
	transform(xyzLigStart);
	RotationTranslation transform_new;
	e_new(i_trans) = refineOne(xyzLigStart,transform_new);
	transforms_new(i_trans) = transform_new*transform;

      }

    }

  }; // class Refinement

} // namespace PRODDL

#endif // AT_PRODDL_REFINE_RIGID_MACOPT_H__
