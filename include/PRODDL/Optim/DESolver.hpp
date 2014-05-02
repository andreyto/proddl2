//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
// Differential Evolution Solver Class
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
// Written By: Lester E. Godwin
//             PushCorp, Inc.
//             Dallas, Texas
//             972-840-0208 x102
//             godwin@pushcorp.com
// Created: 6/8/98
// Last Modified: 6/8/98
// Revision: 1.0

#if !defined(_DESOLVER_H)
#define _DESOLVER_H

#include <string>

namespace PRODDL {

  const int stBest1Exp			= 0;
  const int stRand1Exp			= 1;
  const int stRandToBest1Exp	= 2;
  const int stBest2Exp			= 3;
  const int stRand2Exp			= 4;
  const int stBest1Bin			= 5;
  const int stRand1Bin			= 6;
  const int stRandToBest1Bin	= 7;
  const int stBest2Bin			= 8;
  const int stRand2Bin			= 9;

  class Strategy;
  class StrategyMember;

  class DESolver
  {
  public:

    typedef void (DESolver::*StrategyFunction)(int);

    friend class Strategy;

    DESolver(int dim,int popSize);
    virtual ~DESolver(void);
	
    // Setup() must be called before solve to set min, max, strategy etc.
    void Setup(double min[],double max[], 
	       Strategy* pStrategy,
	       double diffScale=0.7,
	       double crossoverProb=0.9,
	       double scaleFading=0.9);

    void Setup(double min[],double max[], 
	       StrategyFunction pStrategy,
	       double diffScale=0.7,
	       double crossoverProb=0.9,
	       double scaleFading=0.9);


    // Seed the initial population
    virtual void seed();

    // Solve() returns true if EnergyFunction() returns true.
    // Otherwise it runs maxGenerations generations and returns false.
    virtual bool Solve(int maxGenerations);

    // EnergyFunction must be overridden for problem to solve
    // testSolution[] is nDim array for a candidate solution
    // setting bAtSolution = true indicates solution is found
    // and Solve() immediately returns true.
    virtual double EnergyFunction(double testSolution[],bool &bAtSolution) = 0;

    virtual bool checkConstraints(double testSolution[]) { return true; }
	
    int Dimension(void) { return(nDim); }
    int Population(void) { return(nPop); }

    // Call these functions after Solve() to get results.
    double Energy(void) { return(bestEnergy); }
    double *Solution(void) { return(bestSolution); }

    int Generations(void) { return(generations); }

  protected:
    void SelectSamples(int candidate,int *r1,int *r2=0,int *r3=0,
		       int *r4=0,int *r5=0);
    double RandomUniform(double min,double max);

    double& Element(double a[],int b,int c)  { 
      return a[b*nDim+c]; 
    }

    double* RowVector(double a[],int b)  { 
      return &a[b*nDim]; 
    }

    void CopyVector(double a[], double b[]) { 
      for(int i=0; i < nDim; i++) a[i] = b[i]; 
    }


    int nDim;
    int nPop;
    int generations;

    Strategy* m_pStrategy;
    bool m_ownStrategy;

    double scale;
    double probability;
    double scale_rand;

    double trialEnergy;
    double bestEnergy;

    double *trialSolution;
    double *bestSolution;
    double *popEnergy;
    double *population;
    double *min;
    double *max;

    double scaleFading;

  public:
    void Best1Exp(int candidate);
    void Rand1Exp(int candidate);
    void RandToBest1Exp(int candidate);
    void Best2Exp(int candidate);
    void Rand2Exp(int candidate);
    void Best1Bin(int candidate);
    void Rand1Bin(int candidate);
    void RandToBest1Bin(int candidate);
    void Best2Bin(int candidate);
    void Rand2Bin(int candidate);

    void Best1ExpConstr(int candidate);
    void Best2ExpConstr(int candidate);

  private:

    // prevent copy and copy-construction
    DESolver(const DESolver& x) {}
    DESolver& operator= (const DESolver& x) {return *this;}

  }; //class DESolver

  class Strategy {

  public:

    Strategy() {}

    Strategy(DESolver *pDESolver):
      m_pDESolver(pDESolver)
    {}

    void operator() (int candidate) {
      apply(candidate);
    }

    virtual void apply(int candidate) = 0;

  protected:

    DESolver * m_pDESolver;

  }; // class Strategy

  class StrategyMember : public Strategy {

  public:

    StrategyMember() {}

    StrategyMember(DESolver* pDESolver, 
		   DESolver::StrategyFunction strategyFunc):
      Strategy(pDESolver),
      m_strategyFunc(strategyFunc)
    {}

    virtual void apply(int candidate) {
      ((*m_pDESolver).*m_strategyFunc)(candidate);
    }

  protected:

    DESolver::StrategyFunction m_strategyFunc;

  }; // class StrategyMember


    // Combine various parameters and tolerance stop values of DE here in order to
    // provide easy initialization to default values

    struct DEParams {

      int dim;
      int popSize;
      int maxGenerations;
      int maxFuncEvals;
      double cutoffEnergy;
      double deltaEnergy;
      int testGenerations;
      double diffScale;
      double diffScaleRand;
      double crossoverProb;
      double scaleFading;
      std::string strategy;

      DEParams(int dim_):
	dim(dim_),
	popSize(dim*10),
	maxGenerations(popSize*10),
	maxFuncEvals(maxGenerations*popSize),
	cutoffEnergy(-1e60),
	deltaEnergy(1e-7),
	testGenerations(20),
	diffScale(0.7),
	diffScaleRand(diffScale),
	crossoverProb(0.9),
	scaleFading(0.9),
	strategy("Best1ExpConstr")
      {
      }

    };


  class DE : public DESolver {


  public:

    typedef DESolver Base;

    enum {StopUnknown=0,StopLimitGenerations,StopLimitFuncEvals,StopSmallFunc,StopSmallFuncChange,StopEnd};

    DE(const DEParams& params);

    void Setup(double min[],
	       double max[], 
	       StrategyFunction pStrategy = &DESolver::Best1ExpConstr);

    virtual double EnergyFunction(double testSolution[],bool &bAtSolution);

    virtual bool Solve();

    virtual double f(double testSolution[]) = 0;

    int whyStopped() {
      return m_reasonStop;
    }

    const char* whyStoppedStr() {
      int iwhy = whyStopped();
      if( iwhy < 0 || iwhy >= StopEnd ) {
	iwhy = StopUnknown;
      }
      return s_stopStr[iwhy];
    }
      

    DEParams& params() {
      return m_params;
    }

  protected:

    // String values to be indexed by StopXX enums

    static const char* s_stopStr[];

    DEParams m_params;

    int m_testGenerations;
    int m_count;
    int m_reasonStop;
    int m_generation;
    double m_testEnergy;

  }; // class DE

} // namespace PRODDL

#endif // _DESOLVER_H
