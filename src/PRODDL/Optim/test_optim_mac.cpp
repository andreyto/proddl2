//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/to_string.hpp"

#include "PRODDL/Optim/optim_mac.hpp"

#include "PRODDL/Optim/objective_func_bz.hpp"

#include <valarray>

#include <boost/smart_ptr.hpp>

#include "PRODDL/Testing/ctestc.hpp"

CTEST_MODULE()


  // f = sum_i(a*x[i]^2); where a = 1.0/(i+1)
  // min == 0
  void quadratic(const double *p, double *g, int n) {
    for(int i=0; i < n; i++) {
      double a = 1.0/(i+1);
      g[i] = 2*a*p[i];
    }
  }

using namespace blitz;

using namespace PRODDL::Optim;

struct BzPot {


  typedef TinyVector<double,3> Point;
  typedef Array<Point,1> Points;

  typedef ObjectiveFuncContigBzArr<double,Point,BzPot,void> OptimAdaptor;

  void grad(const Points& pos, Points& grad) {
    int n;
    double *_pos = OptimAdaptor::num_cast(pos,n);
    double *_grad = OptimAdaptor::num_cast(grad,n);
    quadratic(_pos,_grad,n);
  }
  
};

int main_test() {
  using namespace PRODDL;
  const int np = 12;
  std::valarray<float> p(np);
  p = 20000;
  //Optim::ObjectiveFuncD1FromPtr<double,void (*)(const double*,double*,int)> f(&quadratic);
  Optim::ObjectiveFuncD1FromPtr<double,void> f(&quadratic);
  Optim::OptimConjGradMac optimizer(f,np);
  // the following does not compile under gcc 2.96:
  //Optim::OptimConjGradMac optimizer(Optim::makeObjectiveFuncD1FromPtr(&quadratic),np);
  boost::scoped_array<float> a(new float[np]);
  for(int i=0; i<np; i++) a[i] = p[i];
  optimizer.optimize(a.get(),np);
  for(int i=0; i<np; i++) p[i] = a[i];
  float r = sqrt((p*=p).sum()/np);
  ATDBGOUT << "RMSD of optimization result from correct minimum is: " << r << std::endl;
  if( r > 1e-6) { throw PRODDL::Testing::test_error("RMSD from true minimum is too large."); }
  BzPot pot;
  BzPot::OptimAdaptor bz_adaptor(pot,&BzPot::grad);
  Optim::OptimConjGradMac optimizer2(bz_adaptor,np);
  p = 20000;
  for(int i=0; i<np; i++) a[i] = p[i];
  optimizer2.optimize(a.get(),np);
  for(int i=0; i<np; i++) p[i] = a[i];
  r = sqrt((p*=p).sum()/np);
  ATDBGOUT << "RMSD of optimization result from correct minimum is: " << r << std::endl;
  if( r > 1e-6) throw PRODDL::Testing::test_error("RMSD from true minimum is too large.");  
  return 0;
}

