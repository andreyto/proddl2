//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/potentialRes2Res.hpp"

#include "PRODDL/exceptions.hpp"

#include "PRODDL/Common/string_util.hpp"

//This potential is basically a
//21x21 matrix with the value of
//potential for each residue-residue
//combination. Columns and rows with C indices from 1 to 20
//correspond to residue types in alphabetical 
//order of three letter names. Row 0 and column 0 correspond
//to residues with some unusual type. Their values are set to
//the average of the corresponding row or column. Cell [0,0]
//is set to the average of all [1:20,1:20] cells.

//Include the matrices Generated from Fabian's data at 6 and 7A

#include "PRODDL/res2res.fabian.6A.cpp"
#include "PRODDL/res2res.fabian.7A.cpp"

namespace PRODDL {

  namespace Res2ResDefault = Res2Res::Fab::A6;

  const int PotentialRes2Res::cResTypes = 20;

  const char* PotentialRes2Res::aResNames3[] = 
    {"ALA", "ARG", "ASN", "ASP", "CYS",
     "GLN", "GLU", "GLY", "HIS", "ILE",
     "LEU", "LYS", "MET", "PHE", "PRO",
     "SER", "THR", "TRP", "TYR", "VAL" };


  //Residues are listed in the order of aResNames3: 1 means residue
  //is hydrophobic; They are: Phe, Ala, Trp, Val, Ile, Ley, Met,
  //Pro

  const int PotentialRes2Res::isHydrophobic[] = 
    {1, 0, 0, 0, 0,
     0, 0, 0, 0, 1,
     1, 0, 1, 1, 1,
     0, 0, 1, 0, 1};

  PotentialRes2Res::PotentialRes2Res() :
     cResNames3(sizeof(aResNames3)/sizeof(aResNames3[0])),
     potMatr(

MatrixType(&Res2ResDefault::a[0][0], blitz::shape(cResTypes+1,cResTypes+1), blitz::duplicateData)),
     potDistance(Res2ResDefault::cutoff)
  {
    // We set 0th row and 0th column to the average of corresponding other columns and rows,
    // and set (0,0)th element to the average of all elements [1:size,1:size].
    // Then we reduce every cell by the value of (0,0).
    // This way, the (i,0) element will reflect the probability of i-th residue to be present at the
    // interface; and the average value of the potential on a uniformly random set of residues will be zero.

    int size = potMatr.rows();
    blitz::Range knownIndexes(1,size-1);
    MatrixType knownMatr(potMatr(knownIndexes,knownIndexes));
    potMatr(knownIndexes,0) = blitz::mean(knownMatr,blitz::tensor::j);
    potMatr(0,knownIndexes) = blitz::mean(knownMatr(blitz::tensor::j,blitz::tensor::i),blitz::tensor::j);
    potMatr(0,0) = blitz::mean(knownMatr);
    potMatr -= potMatr(0,0);
  }

  //return index into Residue-to-Residue potential matrix
  //indices of known residues will be from 1 to cResNames3 inclusive;
  //index 0 will correspond to the unknown residue name, consistent
  //with the structure of the matrix

  int PotentialRes2Res::resName3ToInd(const std::string& sResName3) const
  {
 
    std::string sResName3Up(strupper(sResName3));

    if( sResName3Up.size() != 3 )
      throw potential_error("PotentialRes2Res::resName3ToInd(): 3-letter residue name is expected, instead received: '" +
			    sResName3Up + "'");
   
    for(int i = 0; i < cResNames3; i++) {
      if(sResName3Up == aResNames3[i]) return i+1;
    }
    return 0;
  }

  std::string PotentialRes2Res::resIndToName3(int indRes) const {
    int ind = indRes - 1;
    if( ind >= 0 && ind < cResNames3 )
      return aResNames3[indRes - 1];
    else if( ind == -1 )
      return "";
    else
      throw potential_error("PotentialRes2Res::resIndToName3(): residue index s out of range.");
  }



} // namespace PRODDL

