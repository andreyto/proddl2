//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef RES2RES_POT_H_
#define RES2RES_POT_H_

#include <blitz/array.h>

//residue-residue contact potential

namespace PRODDL {

  class PotentialRes2Res {

  public:

    typedef float T_num;

    typedef blitz::Array<T_num,2> MatrixType;

    PotentialRes2Res();

    //return index into internal Residue-to-Residue potential matrix
    //indices of known residues will be from 1 to cResNames3 inclusive;
    //index 0 will correspond to the unknown residue name, consistent
    //with the structure of the matrix

    int resName3ToInd(const std::string& sResName3) const;


    std::string resIndToName3(int indRes) const;


    // Return a value of potential for two contact residues.
    // The potential is a negative of the contact probability matrix.
    // Thus, the more favorable the contact is, the lower (more negative)
    // value will be returned.
    // Arguments: 'indA' - indexes of two residues in contact (
    // indexes must be obtained by a call to resName3ToInd())
    // If one of the residues is of unknown type (index 0), 
    // then the return value will be the average over all 
    // residues in a column (row) corresponding to another residue.
    // If both residues are of unkonwn type, the average over the entire
    // potential matrix will be returned.


    T_num getValue(int indA, int indB) const {
    
      return - potMatr(indA,indB);

    }


    const MatrixType& getMatrix() const {
      return potMatr;
    }

    // Return a distance at which potential matrix was calculated:
    // Two residues were considered to be in contact if their Cbeta - Cbeta
    // distance was less than getDistanceCutoff().

    T_num getDistanceCutoff() const {
      
      return potDistance;

    }

    int getSize() const {
      return potMatr.rows();
    }

  protected:
    
    static const int cResTypes;

    static const char* aResNames3[];

    const int cResNames3;

    static const int iNameCys = 4; //index of Cys in aResNames3[]

    static const int isHydrophobic[];

    MatrixType potMatr;

    T_num potDistance;

  }; // class PotentialRes2Res

} // namespace PRODDL


    /*

DEPRICATED:

//Find a unary residue based potential, which reflects the probability to
//find a given set of residues at interface on receptor and on ligand
//regardless of which residue is in contact with which. Example: hydrophobic
//residues on receptor and hydrophobic residues on ligand.
//Arguments: each argument is an array of size two - for two structures.
//v[i] - coords of C-beta in i-th structure; ires2name[i][j] - relates an index j
//of v[i] element to the index in potential matrix (index is built from 3-letter
//residue names using ResName3ToInd(); isel[i] - list of indices in v[i], which
//define a set of all CONTACT residues (NOTE the difference with same parameter
//for ResToResPot).
//Return: (Number of hydrophobic residue[receptor] / Number of contact 
//residues[receptor])* 1000 + <same for ligand>.
mfl ResUnaryPot(const vVvect& v, const ivect2& ires2name, const ivect2& isel);
    */


#endif //RES2RES_POT_H_
