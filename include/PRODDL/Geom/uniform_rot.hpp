//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_UNIFORM_ROT_H__
# define AT_UNIFORM_ROT_H__

/* 
   Andrei Tovchigrechko (2001)
   Adapted from Jim Arvo's code from GraphicsGems (see below)
   Changes: 'templatized'; removed dependancy on "GraphicsGems.h"
   and it's datatypes.
*/

#include <exception>
#include <string>

#include "PRODDL/Common/c_array.hpp"

#include <cmath>

/*=========================================================================*
 *  R A N D _ R O T A T I O N      Author: Jim Arvo, 1991                  *
 *                                                                         *
 *  This routine maps three values (x[0], x[1], x[2]) in the range [0,1]   *
 *  into a 3x3 rotation matrix, M.  Uniformly distributed random variables *
 *  x0, x1, and x2 create uniformly distributed random rotation matrices.  *
 *  To create small uniformly distributed "perturbations", supply          *
 *  samples in the following ranges                                        *
 *                                                                         *
 *      x[0] in [ 0, d ]                                                   *
 *      x[1] in [ 0, 1 ]                                                   *
 *      x[2] in [ 0, d ]                                                   *
 *                                                                         *
 * where 0 < d < 1 controls the size of the perturbation.  Any of the      *
 * random variables may be stratified (or "jittered") for a slightly more  *
 * even distribution.                                                      *
 *                                                                         *
 *=========================================================================*/




namespace PRODDL { namespace Geom { namespace Transforms {


 class uniform_rot_error : public std::exception {
 protected:
   std::string m_msg;
 public:
   uniform_rot_error(const std::string& msg) throw (): 
     m_msg(msg) {
   }
   virtual const char* what() const  throw () {
     return m_msg.c_str();
   }
   virtual ~uniform_rot_error() throw () {}
 };

 class uniform_rot_grid_array_err : public uniform_rot_error {
 public:
   uniform_rot_grid_array_err(const std::string& msg)  throw ():
     uniform_rot_error(msg)
   {}
 };


 // second (output) argument must point to T_num[9] and will be filled as 3x3 matrix
 // in row-major order.
  template<typename T_num>  void uniformRotation( const T_num x[3], T_num ares[9])
    {
      for(int i = 0; i < 3; i++) 
	if(x[i] < 0 || x[i] > 1) 
	  throw uniform_rot_error("uniformRotation(): x[] arguments must be within [0,1]");
      
      // cast 'ares' to T_num[][3]
     typename PRODDL::c_array_traits<T_num,2,3>::A_type M =PRODDL::c_array_cast<3>(ares);
      
      T_num PiTimes2 = M_PI*2;
      T_num theta = x[0] * PiTimes2; /* Rotation about the pole (Z).      */
      T_num phi   = x[1] * PiTimes2;       /* For direction of pole deflection. */
      T_num z     = x[2] * 2.0;      /* For magnitude of pole deflection. */

      /* Compute a vector V used for distributing points over the sphere  */
      /* via the reflection I - V Transpose(V).  This formulation of V    */
      /* will guarantee that if x[1] and x[2] are uniformly distributed,  */
      /* the reflected points will be uniform on the sphere.  Note that V */
      /* has length sqrt(2) to eliminate the 2 in the Householder matrix. */

      T_num r  = std::sqrt( z );
      T_num Vx = std::sin( phi ) * r;
      T_num Vy = std::cos( phi ) * r;
      T_num Vz = std::sqrt( 2.0 - z );    

      /* Compute the row vector S = Transpose(V) * R, where R is a simple */
      /* rotation by theta about the z-axis.  No need to compute Sz since */
      /* it's just Vz.                                                    */

      T_num st = std::sin( theta );
      T_num ct = std::cos( theta );
      T_num Sx = Vx * ct - Vy * st;
      T_num Sy = Vx * st + Vy * ct;

      /* Construct the rotation matrix  ( V Transpose(V) - I ) R, which   */
      /* is equivalent to V S - R.                                        */

      M[0][0] = Vx * Sx - ct;
      M[0][1] = Vx * Sy - st;
      M[0][2] = Vx * Vz;

      M[1][0] = Vy * Sx + st;
      M[1][1] = Vy * Sy - ct;
      M[1][2] = Vy * Vz;

      M[2][0] = Vz * Sx;
      M[2][1] = Vz * Sy;
      M[2][2] = 1.0 - z;   /* This equals Vz * Vz - 1.0 */
    }

  // Fill output array 'resM' of T_num[nM*9] elements with uniformly distributed rotation matrices 3x3 produced by
  // varying x[] argument of 'uniformRotation()' in [0,1] interval with step '_step'.
  // int nM - number of elements in output range 'resM'. If nM is zero or negative, this function will not touch
  // the output range, but instead will calculate and return the total number of matrices that would be generated
  // for a given step '_step'. Thus, a typical usage is like this:
  //    nM = uniformRotationGrid(step,0,0);
  //    ... allocate array of nM*9 size ...
  //    uniformRotationGrid(step,aM,nM);

  template<typename T_num> int uniformRotationGrid(T_num _step, T_num* resM, int nM) {

    T_num x[3];
    T_num step[3];

    if(_step < 0) throw uniform_rot_error("uniformRotationGrid(): negative step");

    // later we may decide that step may be different for different coords
    for(int i = 0; i < 3; i++) step[i] = _step;

    int iM = 0;
    for(x[0] = 0; x[0] <= 1; x[0] += step[0])
      for(x[1] = 0; x[1] <= 1; x[1] += step[1])
	for(x[2] = 0; x[2] <= 1; (x[2] += step[2]), iM++) {
	  if(nM>0) {
	    if(iM >= nM) throw uniform_rot_grid_array_err("uniformRotationGrid(): Not enough elements in output array");
	    // cast 'resM' to T_num[][3*3]
	   typename PRODDL::c_array_traits<T_num,2,9>::A_type aM =PRODDL::c_array_cast<9>(resM);
	    uniformRotation(x,aM[iM]);
	  }
	}
    return iM;
  }

} } } // namespace PRODDL { namespace Geom::Transforms

#endif // AT_UNIFORM_ROT_H__

