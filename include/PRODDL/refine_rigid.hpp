//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PRODDL_REFINE_RIGID_H__
#define AT_PRODDL_REFINE_RIGID_H__

// Classes to refines docking predictions by local or semi-local minimization

// Headers for OPT++ library: gradient methods

#include "OPT++/Constraint.h"
#include "OPT++/BoundConstraint.h"
#include "OPT++/OptQNIPS.h"
#include "OPT++/OptCG.h"
#include "OPT++/OptQNewton.h"
#include "OPT++/OptPDS.h"
#include "OPT++/NLF.h"

typedef real OPTPP_real;

// Header for Differential Evolution value-only method

#include "PRODDL/Optim/DESolver.hpp"


//int dummyfunc(real x) { return x; }


#include "PRODDL/potentials.hpp"

#include "PRODDL/Geom/bounding.hpp"

#include "PRODDL/Geom/transformation.hpp"

#include "PRODDL/Geom/move.hpp"

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/g_options.hpp"



namespace PRODDL {

  template<class TargetFunc, typename T_num>
  class OptPPAdaptor : public boost::noncopyable {

  public:

    typedef OptPPAdaptor<TargetFunc,T_num> Self;

    typedef typename Geom::SpaceTraits<T_num>::Point3 Point;
    typedef typename Geom::SpaceTraits<T_num>::VPoint3 Points;

  public:

    OptPPAdaptor(TargetFunc * pTargetFunc,
		 const Points& rbArrCoordsZero, 
		 Points& rbArrCoords):
      m_rbArrGrads(rbArrCoords.copy()),
      m_pTargetFunc(pTargetFunc)
    {
      m_rbArrCoordsZero.reference(rbArrCoordsZero);
      m_rbArrCoords.reference(rbArrCoords);
    }

    T_num optimize() {

      ATLOG_ASSERT_1(m_rbArrCoords.size() == 2);

      const int n_coords = m_rbArrCoords.size()*m_rbArrCoords(0).length();

      const Options& rigOptions = gOptions.getBlock("rigid");

      T_num tolerance;
      rigOptions.getdefault("tolerance",tolerance,1e-6);

      int iterMax;
      rigOptions.getdefault("iterMax",iterMax,200);

      int fevalMax;
      rigOptions.getdefault("fevalMax",fevalMax,iterMax);


      //  OPT++ Nonlinear problem object


      TOLS tol;

      tol.setDefaultTol();

      tol.setStepTol(-1000);
      tol.setFTol(-1000);
      tol.setGTol(tolerance);
      tol.setLSTol(1e-4);
      tol.setMaxStep(0.3);
      tol.setMinStep(1e-15);
      
      tol.setMaxIter(iterMax);
      tol.setMaxFeval(fevalMax);
      tol.setMaxBacktrackIter(iterMax);

      Self * old_self = s_setCurrentSelf(this);

      T_num xyzBound = 0.2; //10 Angstroms
      T_num angBound = 0.1; //25 Degrees 0.44

      Points rbArrLower(2), rbArrUpper(2);
      rbArrLower = m_rbArrCoordsZero;
      rbArrUpper = m_rbArrCoordsZero;
      rbArrLower(0) -= xyzBound;
      rbArrUpper(0) += xyzBound;
      rbArrLower(1) -= angBound;
      rbArrUpper(1) += angBound;

      ColumnVector cvLower(n_coords), cvUpper(n_coords);

      pointsToColumnVector(rbArrLower,cvLower);
      pointsToColumnVector(rbArrUpper,cvUpper);

      //  Create a bound constraint 

      //BoundConstraint bcLower(n_coords, cvLower);
      //BoundConstraint bcUpper(n_coords, cvUpper);

      //CompoundConstraint cc(bcLower,bcUpper);

      Constraint bc = new BoundConstraint(n_coords, cvLower, cvUpper);
      CompoundConstraint cc(bc);

      ////NLF1 nlp(n_coords,s_optGrad,s_optInit,&cc);
      NLF1 nlp(n_coords,s_optGrad,s_optInit);
      ////NLF1 nlp(2,test_rosen,init_test_rosen);
      //NLF0 nlp(n_coords,s_optF,s_optInit);
      //nlp.setDebug();

      ////OptCG optimizer(&nlp,tol);
      OptQNewton optimizer(&nlp,tol);
      ////OptQNIPS optimizer(&nlp,tol);
      ////optimizer.setMeritFcn(ArgaezTapia);

      //optimizer.setSearchStrategy(LineSearch);
      ////PDS:
      //OptPDS optimizer(&nlp,tol);
      //optimizer.setSSS(256);


      //optimizer.setDebug();

      optimizer.setOutputFile("refine_rigid_opt.log",1);

      //optimizer.setIsExpensive(true);

      optimizer.setUpdateModel(s_optUpdateModel);

      optimizer.optimize();

      optimizer.printStatus("RefinementRigid::callOptimizer() finished");

      //optimizer.checkDeriv();

      ColumnVector x_sol = nlp.getXc();
      columnVectorToPoints(x_sol,m_rbArrCoords);      

      T_num f_sol = nlp.getF();

      optimizer.cleanup();

      old_self = s_setCurrentSelf(old_self);

      // Check that no other self has took ownership
      // of the static pointer in the meantime

      ATLOG_ASSERT_1(old_self == this);

      return f_sol;
    }


  protected:

    static Self * s_currentSelf;

    static Self * s_getCurrentSelf() {
      return s_currentSelf;
    }

    static Self * s_setCurrentSelf(Self * self) {
      Self * old_self = s_currentSelf;
      s_currentSelf = self;
      return old_self;
    }

    static
    void s_optInit (int ndim, ColumnVector& x) {
      // set the initial values of x here
      s_getCurrentSelf()->optInit(ndim,x);
    }

    static
    void s_optUpdateModel(int, int,ColumnVector) {
      // can be just empty
    }

    static
    void s_optGrad(int mode, int n, const ColumnVector& x, OPTPP_real& fx, ColumnVector& g, int& result) {
      s_getCurrentSelf()->optGrad(mode,n,x,fx,g,result);
    }

    static
    void s_optF(int n, const ColumnVector& x, OPTPP_real& fx, int& result) {
      s_getCurrentSelf()->optF(n,x,fx,result);
    }


    static
    void pointsToColumnVector(const Points& p, ColumnVector& x) {
      //ColumnVector has unit-based indexes
      for(int i_p = 0, i_x=1; i_p < p.size(); i_p++) {
	for(int j_p = 0; j_p < p(i_p).length(); j_p++, i_x++) {
	  x(i_x) = p(i_p)(j_p);
	}
      }
    }

    static
    void columnVectorToPoints(const ColumnVector& x, Points& p) {
      //ColumnVector has unit-based indexes
      for(int i_p = 0, i_x=1; i_p < p.size(); i_p++) {
	for(int j_p = 0; j_p < p(i_p).length(); j_p++, i_x++) {
	  p(i_p)(j_p) = x(i_x);
	}
      }
    }


    void optInit (int ndim, ColumnVector& x) {
      // set the initial values of x here
      pointsToColumnVector(m_rbArrCoordsZero,x);
    }


    void optGrad(int mode, int n, const ColumnVector& x, OPTPP_real& fx, ColumnVector& g, int& result) {
      // get f, grad or (both?)
      columnVectorToPoints(x,m_rbArrCoords);
      if (mode & NLPFunction) {
	fx  = m_pTargetFunc->f(m_rbArrCoords);
	result = NLPFunction;
      }
      if (mode & NLPGradient) {
	m_pTargetFunc->grad(m_rbArrCoords,m_rbArrGrads);
	pointsToColumnVector(m_rbArrGrads,g);	
	result = NLPGradient;
      }
    }

    void optF(int n, const ColumnVector& x, OPTPP_real& fx, int& result) {
      // get f
      columnVectorToPoints(x,m_rbArrCoords);
      fx  = m_pTargetFunc->f(m_rbArrCoords);
      result = NLPFunction;
    }


    static
    void init_test_rosen (int ndim, ColumnVector& x)
    {
      if (ndim != 2)
	{
	  exit (1);
	}
      x(1) = -1.2;
      x(2) =  1.0;
    }

    static
    void test_rosen(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
    { // Rosenbrock's function
      double f1, f2, x1, x2;

      if (n != 2) return;

      x1 = x(1);
      x2 = x(2);
      f1 = (x2 - x1 * x1);
      f2 = 1. - x1;

      if (mode & NLPFunction) {
	fx  = 100.* f1*f1 + f2*f2;
	result = NLPFunction;
      }
      if (mode & NLPGradient) {
	g(1) = -400.*f1*x1 - 2.*f2;
	g(2) = 200.*f1;
	result = NLPGradient;
      }
    }

  protected:

    Points m_rbArrCoordsZero; // initial coords, should not be changed
    Points m_rbArrCoords; // output coords, reference outside data 2xPoint
    Points m_rbArrGrads; // internal grads, reference own data

    TargetFunc * m_pTargetFunc;

  }; // OptPPAdaptor

  template<class TargetFunc,typename T_num>
  typename OptPPAdaptor<TargetFunc,T_num>::Self * OptPPAdaptor<TargetFunc,T_num>::s_currentSelf = 0;


#define PRODDL_OPT_SELDESTRAT(stratName) \
      else if( strategyName == ""#stratName"" ) { \
	ret = &stratName; \
      }


  template<class TargetFunc, typename T_num>
  class OptDEAdaptor : public DE, public boost::noncopyable {

  public:

    typedef OptDEAdaptor<TargetFunc,T_num> Self;

    typedef typename Geom::SpaceTraits<T_num>::Point3 Point;
    typedef typename Geom::SpaceTraits<T_num>::VPoint3 Points;

    typedef typename common_types::num_vector_type<double>::Type Doubles;

  public:

    OptDEAdaptor(TargetFunc * pTargetFunc,
		 const DEParams& params_):
      DE(params_),
      m_pTargetFunc(pTargetFunc)
    {
      n_coords = params().dim;
      m_dLower.resize(n_coords), m_dUpper.resize(n_coords);
    }


    T_num optimize(const Points& rbArrLower, 
		   const Points& rbArrUpper,
		   Points& rbArrCoords) {

      ATLOG_ASSERT_1(rbArrCoords.size() == 2);

      m_rbArrCoords.reference(rbArrCoords);

      m_rbArrCoordsZero.reference(m_rbArrCoords.copy());

      pointsToDoubles(rbArrLower,m_dLower);
      pointsToDoubles(rbArrUpper,m_dUpper);

      Setup(m_dLower.data(),
	    m_dUpper.data(), 
	    selectStrategy(params().strategy));

      Solve();

      doublesToPoints(Solution(),m_rbArrCoords);      

      T_num f_sol = Energy();

      return f_sol;
    }

    virtual void seed() {

      // This is called from Setup(), which is in turn called from
      // optimized(), so m_dLower, m_dUpper and m_rbArrCoordsZero
      // are already initialized.

      DE::seed();

      Doubles start(nDim);

      pointsToDoubles(m_rbArrCoordsZero,start);

//       for (int i=0; i < nPop; i++) {
// 	bool goodPoint = false;
// 	double e = 0.;
// 	while( ! goodPoint ) {
// 	  for (int j=0; j < nDim; j++) {
// 	    Element(population,i,j) = RandomUniform(min[j],max[j]);
// 	  }
// 	  double e = f(RowVector(population,i));
// 	  if( e < -10 ) {
// 	    goodPoint = true;
// 	  }
// 	}
//       }

      // Replace the first population vector with the initial point

      for(int j = 0; j < nDim; j++) {
	Element(population,0,j) = start(j);
      }

    }

    virtual bool checkConstraints(double testSolution[]) {
      for(int i = 0; i < n_coords; i++) {
	if( testSolution[i] < m_dLower(i) ||
	    testSolution[i] > m_dUpper(i) ) {
	  ATLOG_OUT_4(ATLOGVAR(n_coords) << ATLOGVAR(i) << ATLOGVAR(testSolution[i]) \
		      << ATLOGVAR(m_dLower(i)) << ATLOGVAR(m_dUpper(i)));
	  return false;
	}
      }
      return true;
    }

    virtual double f(double x[]) {
      doublesToPoints(x,m_rbArrCoords);
      return m_pTargetFunc->f(m_rbArrCoords);
    }

  protected:

    DESolver::StrategyFunction selectStrategy(const std::string& strategyName) const {


      struct StSel {
	const char* name;
	DESolver::StrategyFunction p;
	StSel(const char* name_,DESolver::StrategyFunction p_) :
	  name(name_), p(p_) {}
      };

      StSel strats[] = {
	StSel("Best1Exp",&DESolver::Best1Exp),
	StSel("Rand1Exp",&DESolver::Rand1Exp),
	StSel("RandToBest1Exp",&DESolver::RandToBest1Exp),
	StSel("Best2Exp",&DESolver::Best2Exp),
	StSel("Rand2Exp",&DESolver::Rand2Exp),
	StSel("Best1Bin",&DESolver::Best1Bin),
	StSel("Rand1Bin",&DESolver::Rand1Bin),
	StSel("Best1ExpConstr",&DESolver::Best1ExpConstr),
	StSel("Best2ExpConstr",&DESolver::Best2ExpConstr)
      };

      DESolver::StrategyFunction ret = 0;

      for(int i = 0; i < sizeof(strats)/sizeof(strats[0]); i++) {
	if( strats[i].name == strategyName ) {
	  ret = strats[i].p;
	  break;
	}
      }
      
      ATLOG_ASSERT_1(ret != 0);

      return ret;

    }

    static
    void pointsToDoubles(const Points& p, Doubles& x) {
      for(int i_p = 0, i_x=0; i_p < p.size(); i_p++) {
	for(int j_p = 0; j_p < p(i_p).length(); j_p++, i_x++) {
	  x(i_x) = p(i_p)(j_p);
	}
      }
    }

    static
    void doublesToPoints(const Doubles& x, Points& p) {
      for(int i_p = 0, i_x=0; i_p < p.size(); i_p++) {
	for(int j_p = 0; j_p < p(i_p).length(); j_p++, i_x++) {
	  p(i_p)(j_p) = x(i_x);
	}
      }
    }

    static
    void doublesToPoints(const double x[], Points& p) {
      for(int i_p = 0, i_x=0; i_p < p.size(); i_p++) {
	for(int j_p = 0; j_p < p(i_p).length(); j_p++, i_x++) {
	  p(i_p)(j_p) = x[i_x];
	}
      }
    }

  protected:

    int n_coords; // number of coordinates to optimize (6 for xyz+angles)
    Points m_rbArrCoords; // output coords, reference outside data 2xPoint
    Points m_rbArrCoordsZero; // saved initial value of m_rbArrCoords

    // bound constraints
    Doubles m_dLower, m_dUpper;

    TargetFunc * m_pTargetFunc;

  }; // OptDEAdaptor




  template<typename T_num>
  class RefinementRigidValue {

  public:

    typedef RefinementRigidValue<T_num> Self;

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



  protected:

    typedef OptDEAdaptor<Self,T_num> TOptDEAdaptor;

    typedef boost::shared_ptr<TOptDEAdaptor> POptDEAdaptor;


    enum { n_coords = 6 };

    PotTotalNB potTotalNB;



    // buffers for coordinate conversions

    Points xyzCoords, xyzGrads;

    // (starting coords - center of starting coords)
    
    Points xyzCoordsZero;

    // center of starting coords

    PointPair rbCoordsZero;

    Points m_rbArrCoordsZero; // 2xPoint
    Points m_rbArrCoords; // coords, will be referenced and changed by the optimizer

    POptDEAdaptor m_pOptimizer;

    Translation m_iniTransform;

  public:

    RefinementRigidValue() {}

    // NOTE: Array params should be accepted by values, so that we can safely pass
    // here temporary objects returned by view_as_blitz() Python view creation functions
    // and such !!!

    RefinementRigidValue(const Points recPoints, const Points ligPoints, const MolForceParams& mfParams) {

      int runTimeLogLevel;

      gOptions.getdefault("logLevel",runTimeLogLevel,ATLOG_LEVEL_1);

      Logger::setRunTimeLevel(runTimeLogLevel);

      ATLOG_OUT_4(ATLOGVAR(ATLOG_LEVEL) << ATLOGVAR(Logger::getRunTimeLevel()));

      const Options& rigOptions = gOptions.getBlock("rigid");

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

      m_rbArrCoords.resize(2);
      m_rbArrCoordsZero.resize(2);
      m_rbArrCoordsZero = Point(0,0,0);
      xyzCoordsZero.reference(Points(xyzCoords.size()));

      DEParams dePar(n_coords);

      rigOptions.getdefault("tolerance",dePar.deltaEnergy,dePar.deltaEnergy);

      rigOptions.getdefault("iterMax",dePar.maxGenerations,dePar.maxGenerations);

      const Options& deOptions = rigOptions.getBlock("de");

      deOptions.getdefault("nPopulation",dePar.popSize,dePar.popSize);

      deOptions.getdefault("testGenerations",dePar.testGenerations,dePar.testGenerations);

      rigOptions.getdefault("fevalMax",dePar.maxFuncEvals,dePar.popSize*dePar.maxGenerations);

      deOptions.getdefault("diffScale",dePar.diffScale,dePar.diffScale);

      deOptions.getdefault("crossoverProb",dePar.crossoverProb,dePar.crossoverProb);

      deOptions.getdefault("scaleFading",dePar.scaleFading,dePar.scaleFading);

      deOptions.getdefault("strategy",dePar.strategy,dePar.strategy);

      m_pOptimizer.reset(new TOptDEAdaptor(this,dePar));

    }

    T_num f(const Points& rbArrCoords) {

      xyzCoords = xyzCoordsZero;
      // rbArrCoords must contain (xyz,angles) in that order,
      // but RotationTranslation ctor has the order historically reversed
      RotationTranslation tr = m_iniTransform.inverse() * RotationTranslation(rbArrCoords(1),rbArrCoords(0)) * m_iniTransform;
      tr(xyzCoords);
      return potTotalNB.f(xyzCoords);
    }

    T_num 
    refineOne(const Points xyzLigStart,RotationTranslation& resultTransform) {

      const Options& rigOptions = gOptions.getBlock("rigid");

      ATLOG_ASSERT_1(xyzCoordsZero.size() == xyzLigStart.size());      
      xyzCoordsZero = xyzLigStart;

      m_iniTransform = Translation(Point(-1*blitz::mean(xyzCoordsZero)));

      m_rbArrCoords = m_rbArrCoordsZero;

      T_num xyzBound;
      rigOptions.getdefault("xyzBound",xyzBound,1);

      T_num Pi = T_num(4.0)*std::atan(T_num(1.0));
      T_num degToRad = Pi/180.;

      T_num angBound;
      rigOptions.getdefault("angBound",angBound,20);
      angBound *= degToRad;

      Points rbArrLower(2), rbArrUpper(2);
      rbArrLower = m_rbArrCoordsZero;
      rbArrUpper = m_rbArrCoordsZero;
      rbArrLower(0) -= xyzBound;
      rbArrUpper(0) += xyzBound;
      rbArrLower(1) -= angBound;
      rbArrUpper(1) += angBound;

      T_num e = m_pOptimizer->optimize(rbArrLower,rbArrUpper,m_rbArrCoords);

      resultTransform = m_iniTransform.inverse() * RotationTranslation(m_rbArrCoords(1),m_rbArrCoords(0)) * m_iniTransform;

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

  }; // class RefinementRigidValue


  template<typename T_num>
  class RefinementRigidGrad {

  public:

    typedef RefinementRigidGrad<T_num> Self;

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



  protected:

    typedef OptPPAdaptor<Self,T_num> TOptPPAdaptor;

    typedef boost::shared_ptr<TOptPPAdaptor> POptPPAdaptor;

    enum { n_coords = 6 };

    PotTotalNB potTotalNB;



    // buffers for coordinate conversions

    Points xyzCoords, xyzGrads;

    // (starting coords - center of starting coords)
    
    Points xyzCoordsZero;

    // center of starting coords

    PointPair rbCoordsZero;

    Points m_rbArrCoordsZero; // 2xPoint
    Points m_rbArrCoords; // coords, will be referenced and changed by the optimizer
    Points m_rbArrGrads; // output grads, will be referenced and changed by the optimizer

    POptPPAdaptor m_pOptimizer;

  public:

    RefinementRigidGrad() {}

    // NOTE: Array params should be accepted by values, so that we can safely pass
    // here temporary objects returned by view_as_blitz() Python view creation functions
    // and such !!!

    RefinementRigidGrad(const Points recPoints, const Points ligPoints, const MolForceParams& mfParams) {

      int runTimeLogLevel;

      gOptions.getdefault("logLevel",runTimeLogLevel,ATLOG_LEVEL_1);

      Logger::setRunTimeLevel(runTimeLogLevel);

      ATLOG_OUT_4(ATLOGVAR(ATLOG_LEVEL) << ATLOGVAR(Logger::getRunTimeLevel()));

      const Options& rigOptions = gOptions.getBlock("rigid");

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
      m_rbArrCoords.resize(2);
      m_rbArrCoordsZero.resize(2);
      xyzCoordsZero.reference(Points(xyzCoords.size()));

      m_pOptimizer.reset(new TOptPPAdaptor(this,m_rbArrCoordsZero,m_rbArrCoords));

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
#if ATLOG_LEVEL > 8
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
	ATOUTVAR(rbCoordsZero-rbCoords); ATOUTENDL();
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
      PointPair rbCoordsResult;
      rbCoordsResult(0) = rbCoordsZero(0);
      rbCoordsResult(1) = rbCoordsZero(1);
#if ATLOG_LEVEL > 8
      T_num e_xyz_start = 0;
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
	e_xyz_start = e0_1;
	T_num delta_e0 = e0 - e0_1;
	ATOUTVAR(e0); ATOUTVAR(e0_1); ATOUTVAR(delta_e0); ATOUTENDL();
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
      m_rbArrCoords(0) = rbCoordsResult(0);
      m_rbArrCoords(1) = rbCoordsResult(1);
      m_rbArrCoordsZero(0) = rbCoordsZero(0);
      m_rbArrCoordsZero(1) = rbCoordsZero(1);      

      T_num e = m_pOptimizer->optimize();

      rbCoordsResult(0) = m_rbArrCoords(0);
      rbCoordsResult(1) = m_rbArrCoords(1);  
      // make the result to be relative to the initial position
      resultTransform = Geom::Transforms::RigidBody<T_num>::transformation(rbCoordsResult)*
	Geom::Transforms::RigidBody<T_num>::transformation(rbCoordsZero).inverse();
      xyzCoords = xyzLigStart;
      resultTransform(xyzCoords);
#if ATLOG_LEVEL > 8
      {
	m_rbArrCoords(0) = rbCoordsResult(0);
	m_rbArrCoords(1) = rbCoordsResult(1);
	T_num e_rb = f(m_rbArrCoords);
	T_num delta_e_opt_e_rb = e - e_rb;
	T_num e_xyz = potTotalNB.f(xyzCoords);
	T_num delta_e_xyz_e_rb = e_xyz - e_rb;
	T_num e_drop = e_xyz_start - e;
	ATOUTVAR(rbCoordsResult); ATOUTVAR(rbCoordsZero-rbCoordsResult); ATOUTENDL(); 
	ATOUTVAR(e); ATOUTENDL();
	ATOUTVAR(e_rb); ATOUTVAR(e_drop); 
	ATOUTVAR(delta_e_opt_e_rb); ATOUTVAR(delta_e_xyz_e_rb); 
	ATOUTENDL();
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


  }; // class RefinementRigidGrad



} // namespace PRODDL

#endif // AT_PRODDL_REFINE_RIGID_H__
