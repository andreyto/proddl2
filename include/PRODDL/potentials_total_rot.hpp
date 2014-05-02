//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PRODDL_POTENTIALS_H__
#  error POTENTIALS_TOTAL_ROT.HPP MUST BE INCLUDED FROM WITHIN POTENTIALS.HPP
#endif

#ifndef AT_PRODDL_POTENTIALS_TOTAL_ROT_H__
#define AT_PRODDL_POTENTIALS_TOTAL_ROT_H__


// Classes to calculate total potential for the muitiple copy
// optimization of side chain rotamers


  // min of LJ 6-12 is at sigma*ljMinCoeff, and the value is -eps
  //const double ljMinCoeff = 1.122462048309373;

  class PotTotalNonBondedRot {

  public:

    enum ACT_FLAG { ACT_UNDEF = 0x0,
		    ACT_IGN = 0x1, 
		    ACT_INS = 0x2, 
		    ACT_CALC = 0x4,
		    ACT_CALC_INS = ACT_CALC | ACT_INS,
		    ACT_CALC_NO = ~ ACT_CALC,
		    ACT_INS_NO = ~ ACT_INS,
		    ACT_DFL = ACT_CALC_INS };

  protected:


    PartPointsFast m_partPoints;

    SpIter m_iterSphere;

    ForceParAtoms m_fParAtom;

    ForceParNBTable m_fParNBTable;

    Matrix m_aceMatr;

    T_num m_aceCutoff2;

    T_num m_weightEAce;

    PotAce potAce;

    // group-group energy is accumulated here

    Matrix m_fMatr;

    // This data is held by reference,
    // the calling code can change values (but not size)
    // between energy evaluations.
    // All coordinates.

    Points m_points;

    // Same as above regarding change by reference.
    // Group index of every point

    Ints m_indGroups;

    // Same as above regarding change by reference.

    // It should correspond to m_ignMatr, for which
    // values from here serve as indexes.

    // Ignore group index of every point.

    Ints m_indIgn;

    // Same as above regarding change by reference.

    Uints m_indAct;

    // matrix goupIndex x groupIndex to mark
    // pairs of groups for which interactions
    // should not be computed

    IntMatrix m_ignMatr;

  public:
    
    PotTotalNonBondedRot() {}

    PotTotalNonBondedRot(const Points& points,
			 const Ints& indGroups,
			 const Uints& indAct,
			 const Ints& indIgn,
			 const IntMatrix& listIgnoreGroupPairs,
			 int nGroups,
			 int nIgn,
			 const MolForceParams& mfParams,
			 const typename Geom::SpaceTraits<T_num>::Point3Pair& bounds,
			 const Options& options)
    {
      init(points,indGroups,indAct,indIgn,listIgnoreGroupPairs,nGroups,nIgn,mfParams,bounds,options);
    }

    void init(const Points& points,
	      const Ints& indGroups,
	      const Uints& indAct,
	      const Ints& indIgn,
	      const IntMatrix& listIgnoreGroupPairs,
	      int nGroups,
	      int nIgn,
	      const MolForceParams& mfParams,
	      const typename Geom::SpaceTraits<T_num>::Point3Pair& bounds,
	      const Options& options) 
    {

      T_num alpha;

      options.getdefault("alpha",alpha,0.01);

      T_num cutoff;

      options.getdefault("cutoff",cutoff,9.0);

      T_num aceCutoff;

      //options.getdefault("aceCutoff",aceCutoff,6.0);

      aceCutoff = cutoff;

      m_aceCutoff2 = aceCutoff * aceCutoff;

      potAce = PotAce(7.,0.15,0.6,aceCutoff);

      options.getdefault("weightEAce",m_weightEAce,1.);

      m_points.reference(points);

      m_indGroups.reference(indGroups);

      m_indAct.reference(indAct);

      m_indIgn.reference(indIgn);

      m_ignMatr.reference(IntMatrix(nIgn,nIgn));

      setIgnorePairs(listIgnoreGroupPairs);

      m_fMatr.reference(Matrix(nGroups,nGroups));

      //TODO: in order to keep the existing MolForceParams class,
      //we just expect all atoms parameters to be in 'molecule' with 
      //index 0, which is the default in mfParams.paramAtoms()

      m_fParAtom = mfParams.paramAtoms();

      ForceParNBTypes _fParNBTypes(mfParams.paramNBTypes());

      m_fParNBTable = _fParNBTypes.table(alpha,cutoff);

      m_aceMatr.reference(mfParams.aceMatr().copy());

      m_aceMatr *= m_weightEAce;

      m_partPoints.init(bounds(0),bounds(1),cutoff);

      // We also can use a cutoff radius here different from
      // the one used in m_partPoints creation:

      m_iterSphere = SpIter(m_partPoints,PosExtractor(m_points));

      resetF();

    }


    // Clear the result energy matrix

    void resetF() {
      m_fMatr = 0;
    }

    void resetPoints() {
      m_partPoints.zapPoints();
    }

    void resetCycle() {
      resetF();
      resetPoints();
    }


    // Return reference to the matrix with group-group
    // energy values.
    // The calling code is allowed to modify its content -
    // the matrix can be quite large, so we are trying to
    // minimize the number of working copies

    Matrix& getF() {
      return m_fMatr;
    }

    int getNGroups() const {
      return m_fMatr.rows();
    }


    int getNIgn() const {
      return m_ignMatr.rows();
    }

    // By default, interactions of every ignore group with itself will be ignored.
    // It is upto the code that sets up the indIgn assignment to make sure that
    // indIgn corresponds to and least indGroups, if the code does not want
    // the indGroups to interact with themselves.

    void
    setIgnorePairs(const IntMatrix& listIgnoreGroupPairs,bool ignoreSelf=true) {

      m_ignMatr = 0;

      int nPairs = listIgnoreGroupPairs.rows();

      for(int i_pair = 0; i_pair < nPairs; i_pair++) {
	
	int i = listIgnoreGroupPairs(i_pair,0);
	int j = listIgnoreGroupPairs(i_pair,1);

	m_ignMatr(i,j) = 1;
	m_ignMatr(j,i) = 1;

      }

      if( ignoreSelf ) {
	for(int i = 0; i < m_ignMatr.rows(); i++) {
	  m_ignMatr(i,i) = 1;
	}
      }

    }

    void f() {

      //TODO: think how avoid looping over group-group interactions
      //between rotamers of the same residue. So far the only solution
      //seems to be a deep change in the PartPoints object - storing
      //a list of parts' lists instead of a list of points.
      //The generic solution is to use a storage object for each cell,
      //that storage object returning an iterator over objects for
      //a specific input object. Thus, the storage object can do
      //internal memory management and output filtering based on the properties of
      //the elements.
      for(int i_point = 0; i_point < m_points.rows(); i_point++) {
	const Point& point_i = m_points(i_point);
	int i_gr = m_indGroups(i_point);
	unsigned i_act = m_indAct(i_point);
	if( ! (i_act & ACT_IGN) )  { // not ignored 
	  if( i_act & ACT_CALC ) {
	    int i_ign = m_indIgn(i_point);
	    for(m_iterSphere.init(point_i); m_iterSphere.not_end_lazy(); m_iterSphere.next()) {
	      int j_point = *m_iterSphere;
	      int j_gr = m_indGroups(j_point);
	      // We can call this method several times w/o 'reset()'
	      // and disable/enable ACT_CALC in 'm_indAct' in between,
	      // so we still need to check that 'j_act' has flag ACT_CALC here.
	      // ACT_IGN is only checked in the external loop - once the point has
	      // been inserted, it does not matter what its state of ACT_IGN later.
	      unsigned j_act = m_indAct(j_point);
	      if( j_act & ACT_CALC ) {
		int j_ign = m_indIgn(j_point);
		if( ! m_ignMatr(i_ign,j_ign) ) {
		  if( m_iterSphere.check_dist() ) {
		    T_num r2 = m_iterSphere.r2();
		    T_num f = _parNBTableEntry(i_point,j_point).softCoreLJ.f2(r2);
		    if( r2 < m_aceCutoff2 ) {
		      f += aceValue(i_point,j_point) * potAce.f2(r2);
		    }
		    // we do not know in what order i_gr j_gr will appear in the future,
		    // so we need to accumulate in both segments of the matrix
		    m_fMatr(i_gr,j_gr) += f;
		    m_fMatr(j_gr,i_gr) += f;
		  }
		}
	      }
	    }
	  }
	  if( i_act & ACT_INS ) {
	    m_iterSphere.insert(i_point);
	  }
	}
      }
    }

    // Gradient calculation is not implemented yet
    void g(Points& grad);

  protected:

    const ForceParNBTableEntry&
    _parNBTableEntry(int indAtom1,int indAtom2) {
      return m_fParNBTable(m_fParAtom.iType(indAtom1),m_fParAtom.iType(indAtom2));
    }

    T_num aceValue(int indAtom1,int indAtom2) {
      return m_aceMatr(m_fParAtom.aceType(indAtom1),m_fParAtom.aceType(indAtom2));
    }

  };


#endif // AT_PRODDL_POTENTIALS_TOTAL_ROT_H__
