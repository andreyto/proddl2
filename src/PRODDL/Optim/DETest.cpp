//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
// Differential Evolution Test Program
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
// Written By: Lester E. Godwin
//             PushCorp, Inc.
//             Dallas, Texas
//             972-840-0208 x102
//             godwin@pushcorp.com
// Created: 6/8/98
// Last Modified: 6/8/98
// Revision: 1.0

#include <iostream>
#include "PRODDL/Optim/DESolver.hpp"

#include "PRODDL/Common/logger.hpp"

using namespace PRODDL;

// Polynomial fitting problem
class PolynomialSolver : public DESolver
{
public:
	PolynomialSolver(int dim,int pop) : DESolver(dim,pop), count(0) {;}
	double EnergyFunction(double trial[],bool &bAtSolution);

private:
	int count;
};

double PolynomialSolver::EnergyFunction(double *trial,bool &bAtSolution)
{
	int i, j;
	int const M=60;
	double px, x=-1, dx=M, result=0;

	dx = 2.0 / dx;
	for (i=0; i<=M; i++)
	{
		px = trial[0];
		for (j=1;j<nDim;j++)
			px = x*px + trial[j];

		if (px<-1 || px>1)
			result += (1 - px) * (1 - px);

		x += dx;
	}

	px = trial[0];
	for (j=1;j<nDim;j++)
		px = 1.2*px + trial[j];

	px = px - 72.661;
	if (px<0)
		result += px * px;

	px = trial[0];
	for (j=1; j<nDim; j++)
		px = -1.2*px + trial[j];

	px = px - 72.661;
	
	if (px<0)
		result+=px*px;

	if (count++ % nPop == 0)
		printf("%d %lf\n",count / nPop + 1,Energy());
	
	return(result);
}

class RosenSolver : public DE
{
public:
  RosenSolver(const DEParams& params) : DE(params) {}
  double f(double x[]);
};

double RosenSolver::f(double x[])
{
  double y = 0;
  for(int i = 1; i < nDim; i++) {
    double z = (x[i] - x[i-1]*x[i-1]);
    z *= z;
    y += 100. * z;
    z = 1 - x[i-1];
    z *= z;
    y += z;
  }
  return y;
}

namespace PRODDL {

  Logger gLogger;

} // namespace PRODDL


#define N_DIM 9
#define N_POP 100
#define MAX_GENERATIONS	8000

int main()
{

  Logger::setRunTimeLevel(3);

  ATLOG_OUT_1(ATLOGVAR(ATLOG_LEVEL) << ATLOGVAR(Logger::getRunTimeLevel()) <<"\n");

  std::cout << "ATLOG_LEVEL=" << ATLOG_LEVEL << "\t" << Logger::getRunTimeLevel() << "\n";

  double min[N_DIM];
  double max[N_DIM];
  int i;

  DEParams params(N_DIM);

  RosenSolver solver(params);

  for (i=0;i<N_DIM;i++)
    {
      max[i] =  100.0;
      min[i] = -100.0;
    }

  solver.Setup(min,max,
	       //&DESolver::Best2ExpConstr,
	       &DESolver::Best1ExpConstr);

	
  std::cout << "Calculating...\n\n";
  solver.Solve();

  double *solution = solver.Solution();

  std::cout << "\n\nEnergy: " << solver.Energy() << "\n";
  std::cout << "\n\nwhyStopped: " << solver.whyStoppedStr() << "\n";
  std::cout << "\n\nBest Coefficients:\n";
  for (i=0;i<N_DIM;i++)
    std::cout << i << "  :  " << solution[i] << "\n";

  return 0;
}

