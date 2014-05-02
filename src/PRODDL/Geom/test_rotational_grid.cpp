//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/Geom/rotational_grid.hpp"

#include "PRODDL/config.hpp"

#include "PRODDL/Common/debug.hpp"

#include "gtest/gtest.h"


typedef double T_num;

using namespace PRODDL;

using namespace PRODDL::Geom;


///@todo pass parameters for angles file
TEST(RotationalGridTest, Load) {

    std::string dirname;
    
    Parameters().getValue("Params.AngleGridDir",dirname);

    TransformationTraits<T_num>::VRotation rotations = loadRotationalGrid<T_num>(dirname,10);

    for(int i=0; i < rotations.size(); i++)
      std::cout << rotations(i).getTensor() << std::endl;

    int angleFound = 0;

    rotations.reference(loadRotationalGrid<T_num>(dirname,10,angleFound));

    std::cout << "angleFound= " << angleFound << "\n";

    for(int i=0; i < rotations.size(); i++)
      std::cout << rotations(i).getTensor() << std::endl;
    
}
