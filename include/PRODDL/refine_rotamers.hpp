//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PRODDL_REFINE_ROTAMERS_H__
#define AT_PRODDL_REFINE_ROTAMERS_H__


#include <cmath>

#include <boost/utility.hpp>

#include "PRODDL/types.hpp"

#include "PRODDL/potentials.hpp"

#include "PRODDL/Geom/bounding.hpp"

#include "PRODDL/Geom/transformation.hpp"

#include "PRODDL/Geom/move.hpp"

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/g_options.hpp"

namespace PRODDL {


  // Class that implements the computationaly intensive part of the refinement of side chain rotamers 
  // by the multiple copy mean field approach.
  // See the corresponding Python code in Refinement.py for a description of the method and the 
  // literature reference.


  template<typename T_num>
  class RefinementRotamers : boost::noncopyable {

  public:

    typedef typename Potentials<T_num>::PotTotalNonBondedRot Pot;
    typedef typename Potentials<T_num>::MolForceParams MolForceParams;
    typedef typename Geom::SpaceTraits<T_num>::Point3 Point;
    typedef typename Geom::SpaceTraits<T_num>::VPoint3 Points;
    typedef typename Types<T_num>::Matrix Matrix;
    typedef typename Types<T_num>::IntMatrix IntMatrix;
    typedef typename Types<T_num>::Floats Floats;
    typedef typename Types<T_num>::Ints Ints;
    typedef typename Types<T_num>::Uints Uints;
    typedef Geom::Bounding::Box<T_num> BoundingBox;
    typedef Geom::Bounding::Diameter<T_num> BoundingDiameter;
    typedef typename BoundingBox::PointPair PointPair;


  public:


    // The constructor. See definitions of the class members below for a description
    // of each argument (arguments are copied into corresponding class members), except
    // for the 'mfParams' that represents the force field and is passed directly to the
    // object that calculates the total potential.

    // Note:
    // blitz::Array parameters we accept by value, so that we can pass
    // here temporary view objects of Numeric arrays received from Python code,
    // The data inside all parameters is input-only (will not be modified in this code).

    RefinementRotamers(const Points points, 
		       const Ints indGroups, 
		       const Ints indParts,
		       const Ints indMols,
		       const Ints indIgn,
		       const Ints groupPart,
		       const Floats groupFSelf,
		       const IntMatrix rangeMols,
		       const MolForceParams mfParams);





    // Optimize the rotamers.
    //
    // Params:
    // (I  ) points: input coordinates of every possible conformation for all atoms
    // (I  ) isPartVar: flag for variable parts
    // (I/O) groupState: initial/final state of each group (weight)
    // (I/O) partBestGroup: initial/final index of the best (largest weight) group for each part
    // (  O) pointsBest: coordinates of the finally generated best conformation (of the number of real atoms)
    //
    // Note:
    // blitz::Array parameters we accept by value, so that we can pass
    // here temporary view objects of Numeric arrays received from Python code,
    // We can still use methods that change the internal data in output arguments
    // (because the data is held by reference).

    void
    optimize(const Points points,
	     const Ints isPartVar,
	     Floats groupState,
	     Ints partBestGroup,
	     Points pointsBest);



  protected:

    // We hold our own copy of coords - this is safer because it is very easy
    // in the calling Python code to mistakingly reallocate Numpy array

    Points m_points;

    // Object that calculates energy of interaction between groups

    Pot m_potTotal;

    // Constant (independent of mutual receptor-ligand orientation group-group energy)
    // Computed just once

    Matrix m_fConstMatr;

    // group index of every atom (group is something that we want to calculate energy for)
    // This is a work array - its elements are reassigned by this object depending on the stage
    // of calculation

    Ints m_indGroups;

    // number of groups

    int m_nGroups;

    // Action flag for every atom (flags are defined by ACT_XXX enums in potential class)
    // This is a work array - values are reassigned depending on the stage of calculations

    Uints m_indAct;

    // The work copy of m_indAct

    Uints m_indActW;

    // Holds the initial value of m_indGroup passed to the ctor

    Ints m_indGroupsMaster;


    // part index of every atom (part can consist of many groups, an example of a part is 
    // a side chain, that consists of several rotamers)

    Ints m_indParts;

    // molecule index for every atom

    Ints m_indMols;

    // a[group] => part

    Ints m_groupPart;

    // Internal energy of each group.
    // This is precomputed outside amd passed
    // as a parameter to ctor

    Floats m_groupFSelf;

    // indexes of receptor and ligand molecules in
    // m_indMols, and the total number of molecules

    int m_iRec, m_iLig, m_nMols;


    // nAtoms x 2 = [[firstAtom,lastAtom],...] for every molecule

    IntMatrix m_rangeMols;

    // Index of ignore group for each atom
    // If ignPairs[indIgn[atom1],indIgn[atom2]], then the interation between
    // atom1 and atom2 is ignored (e.g. peptide and side chain rotamer)

    Ints m_indIgn;

    // Total number of ignore groups

    int m_nIgn;

    // nAtoms x 2
    // List of pairs of ignore group indexes to ignore during energy calculation,

    IntMatrix m_ignPairs;



  protected:

    Ints flagToIndex(const Ints& flag, int value) {

      Ints ind(blitz::count(flag == value));

      for(int i_ind = 0, i_flag = 0;
	  i_flag < flag.rows();
	  i_flag++) {

	if( flag(i_flag) == value ) {
	  ind(i_ind++) = i_flag;
	}

      }

      return ind;

    }

    void restoreIndGroups() {
      m_indGroups = m_indGroupsMaster;
    }


    // caclulate group interactions that do not depend on the relative positions of
    // receptor and ligand:
    // {R,Rotamers(R)} - Rotamers(R)
    // {L,Rotamers(L)} - Rotamers(L)

    void setupMolecules();


    int countActiveGroups() const {
      
      return blitz::count(!(m_indAct | Pot::ACT_IGN));

    }


    // return the size of array that would be large enough
    // to index with the values of index 'ind',
    // in other words, 'array(ind[i])' would be valid
    // for any 'i', assuming that ind[i] is always >= 0

    int indexSize(const Ints& ind) const {

      return blitz::max(ind) + 1;

    }
    

    void resetPot() {
      m_potTotal.resetCycle();
    }


  }; // class RefinementRotamers


} // namespace PRODDL



  // Implementation


#include "PRODDL/refine_rotamers.cpp"


#endif // AT_PRODDL_REFINE_ROTAMERS_H__
