//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PRODDL_REFINE_ROTAMERS_H__
#  error REFINE_ROTAMERS.CPP MUST BE INCLUDED FROM WITHIN REFINE_ROTAMERS.HPP
#endif

#ifndef AT_PRODDL_REFINE_ROTAMERS_C__
#define AT_PRODDL_REFINE_ROTAMERS_C__



namespace PRODDL {


  // Impelementation of RefineRotamers class

  template<typename T_num>
  RefinementRotamers<T_num>::
  RefinementRotamers(const Points points, 
		     const Ints indGroups, 
		     const Ints indParts,
		     const Ints indMols,
		     const Ints indIgn,
		     const Ints groupPart,
		     const Floats groupFSelf,
		     const IntMatrix rangeMols,
		     const MolForceParams mfParams) {

    // copy original arguments

    m_points.reference(points.copy());
    m_indGroups.reference(indGroups.copy());
    m_indParts.reference(indParts.copy());
    m_indMols.reference(indMols.copy());
    m_indIgn.reference(indIgn.copy());
    m_groupPart.reference(groupPart.copy());
    m_groupFSelf.reference(groupFSelf.copy());
    m_rangeMols.reference(rangeMols.copy());

    // assign derivative member variables

    m_indAct.resize(m_points.rows());
    m_indAct = Pot::ACT_DFL;
    m_indActW.reference(m_indAct.copy());
    m_nGroups = indexSize(m_indGroups);
    m_indGroupsForeman.reference(m_indGroups.copy());

    m_nIgn = indexSize(m_indIgn);
    m_ignPairs.resize(0,2);

    m_iRec = 0;
    m_iLig = 1;
    m_nMols = indexSize(m_indMols);

    int runTimeLogLevel;

    gOptions.getdefault("logLevel",runTimeLogLevel,ATLOG_LEVEL_1);

    Logger::setRunTimeLevel(runTimeLogLevel);

    ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() \
		   << ATLOGVAR(ATLOG_LEVEL) \
		   << ATLOGVAR(Logger::getRunTimeLevel()) <<"\n");


    // select the size of the partiotioning grid as bounding receptor box
    // along the current coordinate axes ('true' as a 2nd argument to
    // BoundingBox ctor)
    // extended in all directions by the diameter of a ligand

    Points recPoints = m_points(blitz::Range(m_rangeMols(m_iRec,0),m_rangeMols(m_iRec,1)));
    Points ligPoints = m_points(blitz::Range(m_rangeMols(m_iLig,0),m_rangeMols(m_iLig,1)));

    BoundingBox boundingBox(recPoints,true);
    PointPair bounds = boundingBox.getDiagonal();
    BoundingDiameter boundingDiameter(ligPoints);
    T_num ligSize = boundingDiameter.getSize();
    bounds[0] -= ligSize;
    bounds[1] += ligSize;

    const Options& rotOptions = gOptions.getBlock("rotamers");

    m_potTotal.init(points,m_indGroups,m_indAct,m_indIgn,m_ignPairs,m_nGroups,m_nIgn,mfParams,bounds,rotOptions);

    setupMolecules();

  }




  template<typename T_num>
  void
  RefinementRotamers<T_num>::
  optimize(const Points points,
	   const Ints isPartVar,
	   Floats groupState,
	   Ints partBestGroup,
	   Points pointsBest) {

    // We first calculate the partMaxGroup 
    // We manipulate m_indGroups (which is referenced inside the m_potTotal)
    // in order to (1) set some parts to constant conformations (all groups but
    // one are inactive for a part), (2) calculate the energy matrix for
    // intermolecular interactions only (intramolecular we precalculated during
    // construction of this object).
    // For that we use groupPart index to see if this group is a preferred state
    // for its part (input value of partBestGroup), and if this group should be
    // treated as a state variable for its part (isPartVar).
    // After all required group-group energies are calculated, we accumulate
    // them into the nVar x nVar energy matrix (where nVar is the number of
    // variable groups, and optimize the weights in partState.
    // We then recalculate the partBestGroup based on the current values
    // in partState
    // Note: Array indStateVar - index of active elements in partState

    // Boltzman law's 'RT' in kJ/mol at T = 300K

    const T_num R = 2.51/300;

    int runTimeLogLevel;

    gOptions.getdefault("logLevel",runTimeLogLevel,ATLOG_LEVEL_1);

    Logger::setRunTimeLevel(runTimeLogLevel);

    const Options& rotOptions = gOptions.getBlock("rotamers");

    T_num T;
    rotOptions.getdefault("T",T,300.);

    const T_num RT = R * T;

    T_num tolerance;
    rotOptions.getdefault("tolerance",tolerance,1e-7);

    T_num lambda;
    rotOptions.getdefault("lambda",lambda,0.5);

    int iterMax;
    rotOptions.getdefault("iterMax",iterMax,3000);


    m_points = points;

    // There are three type of groups: representative states for
    // constant parts, groups of variable parts, and ignored groups

    const int GR_CONST = 1, GR_VAR = 2, GR_IGN = 0;

    Ints grFlag(m_nGroups);

    grFlag = GR_IGN;

    m_indAct = Pot::ACT_CALC_INS;

    // for parts that are not variable, ignore all the alternative groups
    for( int i_pt = 0; i_pt < m_indGroups.rows(); i_pt++) {
      int i_gr = m_indGroups(i_pt);
      int i_part = m_indParts(i_pt);
      int part_best_gr = partBestGroup(i_part);
      int i_mol = m_indMols(i_pt);
      ATLOG_OUT_6(ATLOGVAR(i_gr) << ATLOGVAR(i_part) << ATLOGVAR(part_best_gr) << ATLOGVAR(isPartVar(i_part)));
      if( ! isPartVar(i_part) ) {
	if( i_gr != part_best_gr ) {
	  m_indAct(i_pt) = Pot::ACT_IGN;
	}
	else {
	  grFlag(i_gr) = GR_CONST;
	}
      }
      else {
	grFlag(i_gr) = GR_VAR;
      }
    }

    // indexes of constant and variable groups

    Ints grConst(flagToIndex(grFlag,GR_CONST));

   
    Ints grVar(flagToIndex(grFlag,GR_VAR));

    // indexes of variable parts

    Ints partVar(flagToIndex(isPartVar,1));

    ATLOG_OUT_2("Optimizing " << partVar.rows() << " variable parts having a total of " << grVar.rows() << \
		" states. Number of fixed states is " << grConst.rows());

    // total energies of variable groups

    Floats grVarF(grVar.rows());

    // component of energies of variable groups that
    // comes from interaction with constant groups

    Floats grVarFC(grVar.rows());


    // Zero the energy matrix and empty the PartPoints hash

    resetPot();

    // The following code assumes exactly two molecules (receptor and ligand)

    ATLOG_ASSERT_1(m_nMols == 2);

    m_indActW = m_indAct;

    // Only insert the receptor, do not calculate its interaction with itself,
    // ignore ligand

    m_indAct = blitz::where(m_indMols == m_iRec, m_indAct & Pot::ACT_CALC_NO, m_indAct);
    m_indAct = blitz::where(m_indMols == m_iLig, m_indAct | Pot::ACT_IGN, m_indAct);
    m_potTotal.f();
    m_indAct = m_indActW;

    // Only calculate the interactions for ligand, do not insert (so no interaction with itself),
    // ignore the receptor in this pass (the receptor is already inserted)
    // TODO: we can also ignore the calculation of intermolecular
    // constant-constant part interaction, which we do not need for rotamer optimization -
    // that can be done by splitting the process into several more passes. However, it
    // is probably not that critical, because we are optimizing the interface side chains with
    // a relatively small cutoff value, and therefore the intermolecular constant parts are
    // already removed from each other far enough not to appear in iteration over neighbor cells
    // of a PartPoints grid.

    m_indAct = blitz::where(m_indMols == m_iRec, m_indAct | Pot::ACT_IGN, m_indAct);
    m_indAct = blitz::where(m_indMols == m_iLig, m_indAct | Pot::ACT_INS_NO, m_indAct);
    m_potTotal.f();
    m_indAct = m_indActW;

    // The net result is that only the interaction between receptor and ligand calculated,
    // and no looping is done over receptor-receptor or ligand-ligand positions.

    // Now we have all group-group energies that we need:
    // intramolecular - in m_fConstMatr,
    // intermolecular - in m_potTotal.getF()

    // Because the constant and current energy matrices are complementary, we can just add
    // the constant one to the current one

    Matrix fMatr = m_potTotal.getF();

    fMatr += m_fConstMatr;

    // Steps to accumulate the energies for each rotamer of each variable sidechain (group
    // whose part is isPartVar) and optimize the weights:
    // (1) sum over corresponding row and constant parts columns in fMatr (only one group per
    // constant part), store the sum in a separate vector
    // (2) sum over corresponding row and variable group columns in fMatr with group weights
    // (3) add (1) to (2), recalculate weights, goto (2) till weight self-consistency is reached

    // const component for each variable group, which includes pairwise interaction with
    // all constant groups and internal energy of the variable group itself

    grVarFC = 0;

    for(int i_g = 0; i_g < grVar.rows(); i_g++) {
      int i_grp = grVar(i_g);
      T_num& fc = grVarFC(i_g);
      for(int j_g = 0; j_g < grConst.rows(); j_g++) {
	int j_grp = grConst(j_g);
	fc += fMatr(i_grp,j_grp);
      }
      fc += m_groupFSelf(i_grp);
    }

    Floats groupStateW(groupState.rows());

    Floats partW(isPartVar.rows());


    // variable part and optimization

    int i_iter = 0;

    T_num deltaState = 0;

    // TODO: put grVar elements' centers on the PartPoints object with
    // cutoff larger than max_radius(grVar[i]) + potentialCutoff.
    // Currently we prebuild lists of interacting grVar for each grVar,
    // thus replacing the double loop in all but the first iteration.

    typedef std::vector<int> Vi;
    typedef std::vector<Vi> VVi;

    VVi grVarInter(grVar.rows());

    for(int i = 0; i < grVarInter.size(); i++) {
      grVarInter[i].reserve(100);
    }

    bool madeGrVarInter = false;

    for(i_iter = 0; i_iter < iterMax; i_iter++) {

      grVarF = grVarFC;

      for(int i_g = 0; i_g < grVar.rows(); i_g++) {

	int i_grp = grVar(i_g);

	Vi& inter_i = grVarInter[i_g];

	if( ! madeGrVarInter ) {

	  for(int j_g = 0; j_g < grVar.rows(); j_g++) {
	    int j_grp = grVar(j_g);
	    T_num f = fMatr(i_grp,j_grp);
	    if( f != 0 ) {
	      grVarF(i_g) += groupState(j_grp) * f;
	      inter_i.push_back(j_grp);
	    }
	  }

	}

	else {

	  for(int j_g = 0; j_g < inter_i.size(); j_g++) {
	    int j_grp = inter_i[j_g];
	    grVarF(i_g) += groupState(j_grp) * fMatr(i_grp,j_grp);
	  }

	}
      }

      madeGrVarInter = true;

      groupStateW = 0;

      partW = 0;

      T_num eTot = 0;

      for(int i_g = 0; i_g < grVar.rows(); i_g++) {
	int i_grp = grVar(i_g);
	T_num e = grVarF(i_g);
	T_num w = std::exp( - e / RT );
	groupStateW(i_grp) = w;
	partW(m_groupPart(i_grp)) += w;
	eTot += e;
      }

      T_num S = 0;
      deltaState = 0;

      for(int i_g = 0; i_g < grVar.rows(); i_g++) {
	int i_grp = grVar(i_g);
	T_num g_state = groupState(i_grp);
	T_num g_state_new = lambda * ( groupStateW(i_grp) / partW(m_groupPart(i_grp)) ) +
	  ( 1 - lambda ) * g_state;
	//ATLOG_ASSERT_1(g_state_new <= 1.);
	groupState(i_grp) = g_state_new;
	deltaState += (g_state - g_state_new)*(g_state - g_state_new);
	if( g_state_new > 0 ) {
	  S -= std::log(std::pow(g_state_new,g_state_new));
	}
      }
	
      S /= R;

      deltaState = std::sqrt(deltaState);

      ATLOG_OUT_4("Rotamer refinement: " << ATLOGVAR(i_iter) << ATLOGVAR(deltaState) \
		  << ATLOGVAR(eTot) << ATLOGVAR(S) << "\n");

      if( deltaState < tolerance ) {
	ATLOG_OUT_2("Rotamer matrix tolerance reached at " << ATLOGVAR(i_iter));
	break;
      }

    }

    if( i_iter >= iterMax ) {
      ATLOG_OUT_2("Rotamer matrix iteration limit exceeded at " << ATLOGVAR(deltaState));
    }

    // Find the max state for each part

    // First, set partW to a meaningful state value

    for(int i_g = 0; i_g < grVar.rows(); i_g++) {
      int i_grp = grVar(i_g);
      T_num g_state = groupState(i_grp);
      int i_part = m_groupPart(i_grp);
      partW(i_part) = g_state;
      partBestGroup(i_part) = i_grp;
    }

    // Second, find the max

    for(int i_g = 0; i_g < grVar.rows(); i_g++) {
      int i_grp = grVar(i_g);
      T_num g_state = groupState(i_grp);
      int i_part = m_groupPart(i_grp);
      T_num& max_state = partW(i_part);
      if( g_state > max_state ) {
	max_state = g_state;
	partBestGroup(i_part) = i_grp;
      }
    }


    // Copy points from groups designated as "best representatives"
    // into 'pointsBest'.
    // NOTE: This relies on groups for the same part listing points
    // in the same order, and being the same length

    for(int i_p = 0, i_pb = 0; i_p < points.rows(); i_p++) {
	
      int i_grp = m_indGroups(i_p);

      int i_part = m_indParts(i_p);

      int i_grp_best = partBestGroup(i_part);

      if( i_grp == i_grp_best ) {

	pointsBest(i_pb++) = points(i_p);

      }

    }

    ATLOG_OUT_6(ATLOGVAR(partBestGroup));

  }



  template<typename T_num>
  void
  RefinementRotamers<T_num>::setupMolecules() {

    // We store the calculated values in m_fConstMatr

    m_potTotal.resetCycle();

    // zero matrix of required size

    m_fConstMatr.reference(m_potTotal.getF().copy());

    for( int i_mol=0; i_mol < m_nMols; i_mol++ ) {

      m_indAct = blitz::where(m_indMols == i_mol, Pot::ACT_CALC_INS, Pot::ACT_IGN);

      // forget the previously accumulated data inside m_potTotal

      resetPot();

      m_potTotal.f();

      m_fConstMatr += m_potTotal.getF();

    }

  }

} // namespace PRODDL


#endif // AT_PRODDL_REFINE_ROTAMERS_C__
