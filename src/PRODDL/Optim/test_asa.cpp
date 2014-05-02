//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//


#include "asa_adaptors.h"

#include "atrand.h"

int main()
{
  const int n_atoms = 2000;
  const mfl size_distrib = 20;
  const mfl cutoff = 9.0;
  Vvect av1(n_atoms), av2(n_atoms);
  //generate random points:
  generate(av1.begin(),av1.end(),AtRand::make_gen_ranf(vect(size_distrib/2)));
  generate(av2.begin(),av2.end(),AtRand::make_gen_ranf(vect(size_distrib/2)));
  //object to find energy between teo points:
  vw6_12 fe(3.8,1.0,cutoff,0.1); 
  //optimizer:
  two_body_asa<vw6_12> t(av1,av2,fe,cutoff);
  vect angles(0), xyz(0);
  t.run(angles,xyz,0.5);
  cout << "angles=\t" << angles << endl;
  cout << "xyz=\t" << xyz << endl;
  cout << "cost_function=\t" << t.cost_asa() << endl;
  cout << "exit_code=\t" << t.exit_code_asa() << endl;
  return 0;
} 
