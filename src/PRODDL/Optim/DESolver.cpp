//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/Optim/DESolver.hpp"

#include "PRODDL/Common/logger.hpp"

#include <cmath>

namespace PRODDL {


  DESolver::DESolver(int dim,int popSize) :
    nDim(dim), nPop(popSize),
    generations(0), 
    m_pStrategy(0), m_ownStrategy(false),
    scale(0.7), probability(0.5), 
    scale_rand(0), bestEnergy(0.0),
    trialSolution(0), bestSolution(0),
    popEnergy(0), population(0),
    min(0), max(0),
    scaleFading(0.9)
  {
    trialSolution = new double[nDim];
    bestSolution  = new double[nDim];
    popEnergy	  = new double[nPop];
    population	  = new double[nPop * nDim];
    min           = new double[nDim];
    max           = new double[nDim];
  }

  DESolver::~DESolver(void)
  {
    if (trialSolution) delete trialSolution;
    if (bestSolution) delete bestSolution;
    if (popEnergy) delete popEnergy;
    if (population) delete population;
    if (population) delete population;
    if (min) delete min;
    if (max) delete max;

    trialSolution = bestSolution = popEnergy = population = min = max = 0;

    if( m_ownStrategy ) {
      delete m_pStrategy;
      m_pStrategy = 0;
      m_ownStrategy = false;
    }

  }

  void DESolver::Setup(double *min,double *max,
		       Strategy* pStrategy,
		       double diffScale,
		       double crossoverProb,
		       double scaleFading_)
  {
    int i;

    if( m_ownStrategy ) {
      delete m_pStrategy;
      m_pStrategy = 0;
      m_ownStrategy = false;
    }

    m_pStrategy	= pStrategy;
    scale	= diffScale;
    probability = crossoverProb;
    scaleFading = scaleFading_;

    for(i=0; i < nDim; i++) {
      this->min[i] = min[i];
      this->max[i] = max[i];
    }

    seed();
	
    for (i=0; i < nPop; i++)
      {
	popEnergy[i] = 1.0E20;
      }

    for (i=0; i < nDim; i++) {
      bestSolution[i] = 0.0;
    }

  }

  void DESolver::Setup(double *min,double *max,
		       DESolver::StrategyFunction pStrategy,
		       double diffScale,
		       double crossoverProb,
		       double scaleFading_) {

    Setup(min,max,
	  new StrategyMember(this,pStrategy),
	  diffScale,crossoverProb,scaleFading_);
    m_ownStrategy = true;
  }


  void DESolver::seed() {

    for (int i=0; i < nPop; i++)
      {
	for (int j=0; j < nDim; j++)
	  Element(population,i,j) = RandomUniform(min[j],max[j]);
      }
  }


  bool DESolver::Solve(int maxGenerations)
  {
    int generation;
    int candidate;
    bool bAtSolution;

    bestEnergy = 1.0E20;
    bAtSolution = false;

    for (int i=0; i < nPop; i++)
      {
	double * row = RowVector(population,i);
	double e = EnergyFunction(row,bAtSolution);
	popEnergy[i] = e;
	if( e < bestEnergy ) {
	  bestEnergy = e;
	  CopyVector(bestSolution,row);
	}
      }

    for (generation=0;(generation < maxGenerations) && !bAtSolution;generation++) {
      for (candidate=0; candidate < nPop; candidate++)
	{
	  (*m_pStrategy)(candidate);
	  trialEnergy = EnergyFunction(trialSolution,bAtSolution);

	  if (trialEnergy < popEnergy[candidate])
	    {
	      // New low for this candidate
	      popEnergy[candidate] = trialEnergy;
	      CopyVector(RowVector(population,candidate),trialSolution);

	      // Check if all-time low
	      if (trialEnergy < bestEnergy)
		{
		  bestEnergy = trialEnergy;
		  CopyVector(bestSolution,trialSolution);
		}
	    }
	}
      ATLOG_OUT_3(ATLOGVAR(generation) << ATLOGVAR(trialEnergy) << \
		  ATLOGVAR(bestEnergy) << ATLOGVAR(bAtSolution)<< "\n");
    }

    generations = generation;

    return(bAtSolution);
  }

  void DESolver::Best1Exp(int candidate)
  {
    int r1, r2;
    int n;

    SelectSamples(candidate,&r1,&r2);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
      {
	trialSolution[n] = bestSolution[n]
	  + scale * (Element(population,r1,n)
		     - Element(population,r2,n));

	n = (n + 1) % nDim;
      }

    return;
  }

  void DESolver::Rand1Exp(int candidate)
  {
    int r1, r2, r3;
    int n;

    SelectSamples(candidate,&r1,&r2,&r3);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
      {
	trialSolution[n] = Element(population,r1,n)
	  + scale * (Element(population,r2,n)
		     - Element(population,r3,n))
	  + scale_rand * RandomUniform(min[n],max[n]);
	n = (n + 1) % nDim;
      }

    return;
  }

  void DESolver::RandToBest1Exp(int candidate)
  {
    int r1, r2;
    int n;

    SelectSamples(candidate,&r1,&r2);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
      {
	trialSolution[n] += scale * (bestSolution[n] - trialSolution[n])
	  + scale * (Element(population,r1,n)
		     - Element(population,r2,n));
	n = (n + 1) % nDim;
      }

    return;
  }

  void DESolver::Best2Exp(int candidate)
  {
    int r1, r2, r3, r4;
    int n;

    SelectSamples(candidate,&r1,&r2,&r3,&r4);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
      {
	trialSolution[n] = bestSolution[n] +
	  scale * (Element(population,r1,n)
		   + Element(population,r2,n)
		   - Element(population,r3,n)
		   - Element(population,r4,n));
	n = (n + 1) % nDim;
      }

    return;
  }

  void DESolver::Rand2Exp(int candidate)
  {
    int r1, r2, r3, r4, r5;
    int n;

    SelectSamples(candidate,&r1,&r2,&r3,&r4,&r5);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
      {
	trialSolution[n] = Element(population,r1,n)
	  + scale * (Element(population,r2,n)
		     + Element(population,r3,n)
		     - Element(population,r4,n)
		     - Element(population,r5,n));
	n = (n + 1) % nDim;
      }

    return;
  }

  void DESolver::Best1Bin(int candidate)
  {
    int r1, r2;
    int n;

    SelectSamples(candidate,&r1,&r2);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; i < nDim; i++) 
      {
	if ((RandomUniform(0.0,1.0) < probability) || (i == (nDim - 1)))
	  trialSolution[n] = bestSolution[n]
	    + scale * (Element(population,r1,n)
		       - Element(population,r2,n));
	n = (n + 1) % nDim;
      }

    return;
  }

  void DESolver::Rand1Bin(int candidate)
  {
    int r1, r2, r3;
    int n;

    SelectSamples(candidate,&r1,&r2,&r3);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; i < nDim; i++) 
      {
	if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
	  trialSolution[n] = Element(population,r1,n)
	    + scale * (Element(population,r2,n)
		       - Element(population,r3,n));
	n = (n + 1) % nDim;
      }

    return;
  }

  void DESolver::RandToBest1Bin(int candidate)
  {
    int r1, r2;
    int n;

    SelectSamples(candidate,&r1,&r2);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; i < nDim; i++) 
      {
	if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
	  trialSolution[n] += scale * (bestSolution[n] - trialSolution[n])
	    + scale * (Element(population,r1,n)
		       - Element(population,r2,n));
	n = (n + 1) % nDim;
      }

    return;
  }

  void DESolver::Best2Bin(int candidate)
  {
    int r1, r2, r3, r4;
    int n;

    SelectSamples(candidate,&r1,&r2,&r3,&r4);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; i < nDim; i++) 
      {
	if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
	  trialSolution[n] = bestSolution[n]
	    + scale * (Element(population,r1,n)
		       + Element(population,r2,n)
		       - Element(population,r3,n)
		       - Element(population,r4,n));
	n = (n + 1) % nDim;
      }

    return;
  }

  void DESolver::Rand2Bin(int candidate)
  {
    int r1, r2, r3, r4, r5;
    int n;

    SelectSamples(candidate,&r1,&r2,&r3,&r4,&r5);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; i < nDim; i++) 
      {
	if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
	  trialSolution[n] = Element(population,r1,n)
	    + scale * (Element(population,r2,n)
		       + Element(population,r3,n)
		       - Element(population,r4,n)
		       - Element(population,r5,n));
	n = (n + 1) % nDim;
      }

    return;
  }


  void DESolver::Best1ExpConstr(int candidate)
  {
    int r1, r2;
    int n;

    SelectSamples(candidate,&r1,&r2);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
      {

	double x = trialSolution[n];
	double diff = Element(population,r1,n) - Element(population,r2,n)
	  + scale_rand/scale * RandomUniform(min[n],max[n]);

	double scale_here = scale;
	for( bool ok_constr = false; ! ok_constr; ) {
	  double y = bestSolution[n] + scale_here * diff;
	  trialSolution[n] = y;
	  ATLOG_OUT_4(ATLOGVAR(y) << ATLOGVAR(scale_here) << ATLOGVAR(scaleFading) << "\n");
	  ok_constr = checkConstraints(trialSolution);
	  if( ! ok_constr ) {
	    trialSolution[n] = x;
	    scale_here *= scaleFading;
	    if( scale_here < 1e-30 ) {
	      scale_here = 0.;
	    }
	  }
	}
	n = (n + 1) % nDim;
      }

  }


  void DESolver::Best2ExpConstr(int candidate)
  {
    int r1, r2, r3, r4;
    int n;

    SelectSamples(candidate,&r1,&r2,&r3,&r4);
    n = (int)RandomUniform(0.0,(double)nDim);

    CopyVector(trialSolution,RowVector(population,candidate));
    for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
      {

	double x = trialSolution[n];
	double diff = (Element(population,r1,n)
		       + Element(population,r2,n)
		       - Element(population,r3,n)
		       - Element(population,r4,n));
	double scale_here = scale;
	for( bool ok_constr = false; (! ok_constr) && scale_here >= 1e-10; ) {
	  double y = bestSolution[n] + scale_here * diff;
	  trialSolution[n] = y;
	  ATLOG_OUT_4(ATLOGVAR(y) << ATLOGVAR(scale_here) << ATLOGVAR(scaleFading) << "\n");
	  ok_constr = checkConstraints(trialSolution);
	  if( ! ok_constr ) {
	    trialSolution[n] = x;
	    scale_here *= scaleFading;
	    if( scale_here < 1e-30 ) {
	      scale_here = 0.;
	    }
	  }
	}

	n = (n + 1) % nDim;
      }

  }



  void DESolver::SelectSamples(int candidate,int *r1,int *r2,
			       int *r3,int *r4,int *r5)
  {
    if (r1)
      {
	do
	  {
	    *r1 = (int)RandomUniform(0.0,(double)nPop);
	  }
	while (*r1 == candidate);
      }

    if (r2)
      {
	do
	  {
	    *r2 = (int)RandomUniform(0.0,(double)nPop);
	  }
	while ((*r2 == candidate) || (*r2 == *r1));
      }

    if (r3)
      {
	do
	  {
	    *r3 = (int)RandomUniform(0.0,(double)nPop);
	  }
	while ((*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1));
      }

    if (r4)
      {
	do
	  {
	    *r4 = (int)RandomUniform(0.0,(double)nPop);
	  }
	while ((*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) || (*r4 == *r1));
      }

    if (r5)
      {
	do
	  {
	    *r5 = (int)RandomUniform(0.0,(double)nPop);
	  }
	while ((*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3)
	       || (*r5 == *r2) || (*r5 == *r1));
      }

    return;
  }

/*------Constants for RandomUniform()-----*/
const long SEED = 3;
const long IM1 = 2147483563;
const long IM2 = 2147483399;
const double AM = (1.0/IM1);
const long IMM1 = (IM1-1);
const long IA1 = 40014;
const long IA2 = 40692;
const long IQ1 = 53668;
const long IQ2 = 52774;
const long IR1 = 12211;
const long IR2 = 3791;
const int  NTAB = 32;
const long NDIV = (1+IMM1/NTAB);
const double EPS = 1.2e-7;
const double RNMX = (1.0-EPS);

  double DESolver::RandomUniform(double minValue,double maxValue)
  {
    long j;
    long k;
    static long idum;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double result;

    if (iy == 0)
      idum = SEED;

    if (idum <= 0)
      {
	if (-idum < 1)
	  idum = 1;
	else
	  idum = -idum;

	idum2 = idum;

	for (j=NTAB+7; j>=0; j--)
	  {
	    k = idum / IQ1;
	    idum = IA1 * (idum - k*IQ1) - k*IR1;
	    if (idum < 0) idum += IM1;
	    if (j < NTAB) iv[j] = idum;
	  }

	iy = iv[0];
      }

    k = idum / IQ1;
    idum = IA1 * (idum - k*IQ1) - k*IR1;

    if (idum < 0)
      idum += IM1;

    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;

    if (idum2 < 0)
      idum2 += IM2;

    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = idum;

    if (iy < 1)
      iy += IMM1;

    result = AM * iy;

    if (result > RNMX)
      result = RNMX;

    result = minValue + result * (maxValue - minValue);
    return(result);
  }

  DE::DE(const DEParams& params_):
    DE::Base(params_.dim,params_.popSize),
    m_params(params_),
    m_count(0),
    m_reasonStop(DE::StopLimitGenerations),
    m_generation(0),
    m_testEnergy(0)
  {}

  void DE::Setup(double min[],
		 double max[], 
		 StrategyFunction pStrategy) {
    Base::Setup(min,
		max,
		pStrategy,
		m_params.crossoverProb,
		m_params.scaleFading);
    scale_rand = m_params.diffScaleRand;
    m_count = 0;
    m_reasonStop = StopLimitGenerations;
    m_generation = 0;
    m_testEnergy = 0;
  }

  double DE::EnergyFunction(double testSolution[],bool &bAtSolution) {

    double e = this->f(testSolution);
        
    m_count += 1;

    if( m_count >= m_params.maxFuncEvals ) {
      bAtSolution = true;
      m_reasonStop = StopLimitFuncEvals;
    }
        
    // m_count is per evaluation, m_count % nPop is per this->generation
        
    if( (m_count-1) % nPop == 0) {
      m_generation = m_count / nPop;
      //print self.count, nPop, self.generation, self.e, self.Solution()

      // we will be "done" if the best energy is less than or equal to the
      // cutoff energy (defaults to some very negative number)
            
      if( bestEnergy <= m_params.cutoffEnergy ) {
	bAtSolution = true;
	m_reasonStop = StopSmallFunc;
      }

      // we will be "done" if the best energy is changed by less that deltaEnergy
      // every "testGenerations" generations
            
      if( m_generation == m_params.testGenerations ) {  // set initial test energy
	m_testEnergy = bestEnergy;
      }

      // test every self.testGenerations generations after the initialization above
            
      if( m_generation > m_params.testGenerations &&
	  m_generation % m_params.testGenerations == 0 ) {

	// if energy changes by less than deltaEnergy in "testGenerations"
	// generations, stop
                
	double deltaEnergy = std::fabs(m_testEnergy - bestEnergy);

	double absEnergy = std::fabs(bestEnergy);

	double e_min = 1e60, e_max = -1e60;

	for(int i = 0; i < nPop; i++ ) {

	  e = popEnergy[i];

	  if(e < e_min) {
	    e_min = e;
	  }

	  if(e > e_max) {
	    e_max = e;
	  }

	}

	if(std::fabs(e_max - e_min) < m_params.deltaEnergy) {

	  if( ( absEnergy >  m_params.deltaEnergy && deltaEnergy/absEnergy < m_params.deltaEnergy ) ||
	      ( absEnergy <= m_params.deltaEnergy && deltaEnergy < m_params.deltaEnergy ) ) {
	    //if( scale_rand > 0 ) {
	    //  scale_rand = 0;
	    //}
	    //else {
	    bAtSolution = true;
	    m_reasonStop = StopSmallFuncChange;
	    //seed();
	    //}
	  }

	}
	
	m_testEnergy = bestEnergy;
	
      }

    }
    
    return e;
  }



  bool DE::Solve() {

    return Base::Solve(m_params.maxGenerations);

  }

  // Values must correspond to enums with the same names and indexes from the class definition except StopEnd which is the size of the array
  const char* DE::s_stopStr[StopEnd] = {"StopUnknown","StopLimitGenerations","StopLimitFuncEvals","StopSmallFunc","StopSmallFuncChange"};


} // namespace PRODDL
