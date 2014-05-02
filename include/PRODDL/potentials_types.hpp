//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PRODDL_POTENTIALS_H__
#  error POTENTIALS_TYPES.HPP MUST BE INCLUDED FROM WITHIN POTENTIALS.HPP
#endif

#ifndef AT_PRODDL_POTENTIALS_TYPES_H__
#define AT_PRODDL_POTENTIALS_TYPES_H__


// Classes to build tables of non-bonded pairwise potentials

typedef int IndAtomType;

class ForceParAtoms {

public:


  Points m_pos;

  Floats m_mass;

  Ints m_iType;

  Ints m_aceType;

  ForceParAtoms& operator= (const ForceParAtoms& other) {

    if( this != &other ) {

      m_pos.reference(other.m_pos);

      m_iType.reference(other.m_iType);

      m_aceType.reference(other.m_aceType);

      m_mass.reference(other.m_mass);

    }

    return *this;

  }

  Ints iType() const {

    return m_iType;

  }

  Floats mass() const {

    return m_mass;

  }

  Points pos() const {

    return m_pos;

  }

  int iType(int iAt) const {

    return m_iType(iAt);

  }

  int aceType(int iAt) const {

    return m_aceType(iAt);

  }

  Float mass(int iAt) const {

    return m_mass(iAt);

  }

  Point pos(int iAt) const {

    return m_pos(iAt);

  }


};


class ForceParNBTableEntry {

public:

  SoftCoreLJ softCoreLJ;

};


typedef typename common_types::num_matrix_type<ForceParNBTableEntry>::Type ForceParNBTableArr;


class ForceParNBTable {

protected:

  ForceParNBTableArr table;

public:

  ForceParNBTable() {}

  ForceParNBTable(int _size) :
    table(_size,_size) 
  {}


  ForceParNBTableEntry& operator() (IndAtomType indType1, IndAtomType indType2) {

    return table(indType1,indType2);

  }

  const
  ForceParNBTableEntry& operator() (IndAtomType indType1, IndAtomType indType2) const {

    return table(indType1,indType2);

  }


  ForceParNBTable& operator= (const ForceParNBTable& other) {

    if( this != &other ) {

      table.reference(other.table);

    }

    return *this;

  }

  int size() const {

    return table.size();

  }  
  
};


class ForceParNBTypes {

public:

  Floats sigma;

  Floats eps;

  int mix;

  ForceParNBTypes& operator= (const ForceParNBTypes& other) {

    if( this != &other ) {

      sigma.reference(other.sigma);

      eps.reference(other.eps);

      mix = mix;

    }

    return *this;

  }

  ForceParNBTable table(T_num alpha, T_num cutoff) const {

    int n = size();

    ForceParNBTable _table(n);

    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++) {
	T_num _sigma;
	if( mix == LJ_MIX_0 ) {
	  _sigma = (sigma(i) + sigma(j))/2.;
	}
	else {
	  _sigma = std::sqrt(sigma(i)*sigma(j));
	}
	T_num _eps = std::sqrt(eps(i)*eps(j));
	_table(i,j).softCoreLJ = SoftCoreLJ(_sigma,_eps,alpha,cutoff);
      }
	
    return _table;

  }

  int size() const {

    return sigma.size();

  }

};


class MolForceParams {

public:

  typedef std::vector<ForceParAtoms> ArrForceParAtoms;

  ArrForceParAtoms fpAtoms;

  ForceParNBTypes nbTypes;

  Matrix m_aceMatr;

  ForceParAtoms paramAtoms(int iMol=0) const {

    return fpAtoms[iMol];

  }

  ForceParNBTypes paramNBTypes() const {

    return nbTypes;

  }

  Matrix aceMatr() const {
    return m_aceMatr;
  }

};


#endif // AT_PRODDL_POTENTIALS_TYPES_H__
