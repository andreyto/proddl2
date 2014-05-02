//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT_PRODDL_POTENTIALS_H__
#  error POTENTIALS_TOTAL.HPP MUST BE INCLUDED FROM WITHIN POTENTIALS.HPP
#endif

#ifndef AT_PRODDL_POTENTIALS_TOTAL_H__
#define AT_PRODDL_POTENTIALS_TOTAL_H__


// Classes to calculate total potentials of a system


// min of LJ 6-12 is at sigma*ljMinCoeff, and the value is -eps
//const double ljMinCoeff = 1.122462048309373;


class PotTotalNonBonded {


protected:

	// if !0 - do use harmonic basin
	int harmonicBasinUse;

	Harmonic harmonicBasin;

	PartPoints partPoints;

	Points points1;

	ForceParAtoms fParAtom1, fParAtom2;

	ForceParNBTable fParNBTable;

	ForceParNBTypes m_fParNBTypes;

	Point points1Center;

	//DEBUG:
	SoftCoreLJ softCoreLJ;

	T_num m_cutoffP2;

	T_num m_alpha;

	Matrix m_aceMatr;

	T_num m_aceCutoff2;

	T_num m_weightEAce;

	PotAce potAce;

	Options m_opts;

public:

	PotTotalNonBonded() {}

	PotTotalNonBonded(const Points& _points1,
		const MolForceParams& mfParams,
		const typename Geom::SpaceTraits<T_num>::Point3Pair& bounds,
		const Options& options)
	{
		init(mfParams,bounds,options);
	}

	void init(const Points& _points1,
		const MolForceParams& mfParams,
		const typename Geom::SpaceTraits<T_num>::Point3Pair& bounds,
		const Options& options) 
	{

		m_opts = options;

		//const T_num sigma = 0.4; //0.34;
		//const T_num eps = 0.5; // 0.36;
		//alpha = 0.35;

		T_num alpha;
		options.getdefault("alpha",alpha,0.0);

		m_alpha = alpha;

		// make harmonic basin constant small enough
		T_num harmonicBasinK;
		options.getdefault("harmonicBasinK",harmonicBasinK,1.e-9);
		const T_num k = harmonicBasinK / Math::dotSelf(bounds(1) - bounds(0));

		harmonicBasin = Harmonic(k);

		options.getdefault("harmonicBasinUse",harmonicBasinUse,0);

		points1.reference(_points1.copy());

		fParAtom1 = mfParams.paramAtoms(0);

		fParAtom2 = mfParams.paramAtoms(1);


		points1Center = blitz::mean(points1);

		T_num cutoff;

		options.getdefault("cutoff",cutoff,9.0);

		ATLOG_OUT_4(ATLOGVAR(cutoff));

		m_cutoffP2 = cutoff * cutoff;

		T_num aceCutoff;

		//options.getdefault("aceCutoff",aceCutoff,6.0);

		aceCutoff = cutoff;

		m_aceCutoff2 = aceCutoff * aceCutoff;

		potAce = PotAce(7.,0.15,0.6,aceCutoff);
		//potAce = PotAce(7.,0.15,0.7,aceCutoff);

		options.getdefault("weightEAce",m_weightEAce,1.);

		m_fParNBTypes = mfParams.paramNBTypes();

		fParNBTable = m_fParNBTypes.table(alpha,cutoff);

		m_aceMatr.reference(mfParams.aceMatr().copy());

		m_aceMatr *= m_weightEAce;

		partPoints.init(bounds(0),bounds(1),cutoff);

		partPoints.insert(points1);
	}


	void setPoints1(const Points& _points1) {

		ATLOG_ASSERT_1(points1.rows() == _points1.rows());

		points1 = _points1;

		points1Center = blitz::mean(points1);

		partPoints.zapPoints();

		partPoints.insert(points1);

	}

	T_num f(const Points& points2) {

		VIPair indexPairs;
		fvect distanceP2;

		partPoints.search(points2,indexPairs,distanceP2);

		int n_pairs = distanceP2.size();

		T_num f_total = 0;

#ifdef PRODDL_CONTRIB
		int doCustomEnergy;
		m_opts.getdefault("doCustomEnergy",doCustomEnergy,0);
		if( doCustomEnergy ) {
			T_num delta;
			m_opts.getBlock("customEnergy").get("delta",delta);
			for(int i_pair=0; i_pair < n_pairs; i_pair++) {
				const IPair& ind_pair = indexPairs[i_pair];
				T_num rP2 = distanceP2[i_pair];
				int iType1 = fParAtom1.iType(ind_pair(1));
				int iType2 = fParAtom2.iType(ind_pair(0));
				f_total += DKReal(customEnergy(DKReal(m_alpha), 
					DKReal((m_fParNBTypes.sigma(iType1)+
					m_fParNBTypes.sigma(iType2)) / 2), 
					DKReal(std::sqrt(m_fParNBTypes.eps(iType1)*
					m_fParNBTypes.eps(iType2))), 
					DKReal(delta), 
					DKReal(std::sqrt(rP2))));
			}
			return f_total;
		}
#endif

		for(int i_pair=0; i_pair < n_pairs; i_pair++) {
			const IPair& ind_pair = indexPairs[i_pair];
			T_num rP2 = distanceP2[i_pair];



			//ATTENTION: ind_pair(0) is for points2, ind_pair(1) -- points1
			//DEBUG:
			f_total += _parNBTableEntry(ind_pair(1),ind_pair(0)).softCoreLJ.f2(rP2);
			//f_total += softCoreLJ.f2(rP2);

			if( rP2 < m_aceCutoff2 ) {
				f_total += aceValue(ind_pair(1),ind_pair(0)) * potAce.f2(rP2);
			}
		}

		if( harmonicBasinUse ) {
			for(int i_point = 0; i_point < points2.size(); i_point++) {
				f_total += harmonicBasin.f2(Math::dotSelf(points2(i_point) - points1Center));
			}
		}

		bool testMode = false;

		if( testMode ) {
			T_num f_totalDoubleLoop = fDoubleLoop(points2);
			ATLOG_OUT_1("In test mode: " << ATLOGVAR(f_total) << ATLOGVAR(f_totalDoubleLoop) \
				<< ATLOGVAR(f_total-f_totalDoubleLoop));
			ATLOG_ASSERT_1(Math::are_close(f_total,f_totalDoubleLoop));
			return f_totalDoubleLoop;
		}

		return f_total;

	}


	void g(const Points& points2,Points& grad2) {

		VIPair indexPairs;
		fvect distanceP2;

		partPoints.search(points2,indexPairs,distanceP2);

		int n_pairs = distanceP2.size();

		grad2 = 0;


#ifdef PRODDL_CONTRIB
		int doCustomEnergy;
		m_opts.getdefault("doCustomEnergy",doCustomEnergy,0);
		if( doCustomEnergy ) {
			throw not_supported_error("We do not expect CustomEnergy to implement gradient calculation.");
		}
#endif

		for(int i_pair=0; i_pair < n_pairs; i_pair++) {
			const IPair& ind_pair = indexPairs[i_pair];
			T_num rP2 = distanceP2[i_pair];
			//ATTENTION: ind_pair(0) is for points2, ind_pair(1) -- points1
			//DEBUG:
			Point p = points2(ind_pair(0)) - points1(ind_pair(1));
			grad2(ind_pair(0)) += _parNBTableEntry(ind_pair(1),ind_pair(0)).
				softCoreLJ.g2(rP2,p);
			//grad2(ind_pair(0)) += softCoreLJ.g2(rP2,Point(points2(ind_pair(0)) - points1(ind_pair(1))));
			if( rP2 < m_aceCutoff2 ) {
				grad2(ind_pair(0)) += aceValue(ind_pair(1),ind_pair(0)) * potAce.g2(rP2,p);
			}
		}

		if( harmonicBasinUse ) {
			for(int i_point = 0; i_point < points2.size(); i_point++) {
				grad2(i_point) += harmonicBasin.g2(points2(i_point) - points1Center);
			}
		}

		bool testMode = false;

		if( testMode ) {
			Points grad2DoubleLoop(grad2.size());
			gDoubleLoop(points2,grad2DoubleLoop);
			Point s = Math::dotSelf(grad2DoubleLoop - grad2);
			T_num rmsdGradDoubleLoop = std::sqrt(blitz::sum(s))/T_num(grad2.size());
			T_num normGrad = std::sqrt(blitz::sum(Math::dotSelf(grad2)));
			ATLOG_OUT_1("In test mode: " << ATLOGVAR(normGrad) << ATLOGVAR(rmsdGradDoubleLoop));
			ATLOG_ASSERT_1(rmsdGradDoubleLoop < 1e-7);
			grad2 = grad2DoubleLoop;
		}

	}

	T_num fDoubleLoop(const Points& points2) {


		T_num f_total = 0;

		for(int i = 0; i < points1.size(); i++) {
			for(int j = 0; j < points2.size(); j++) {
				T_num rP2 = Math::dotSelf(points2(j) - points1(i));
				if( rP2 < m_cutoffP2 ) {
					f_total += _parNBTableEntry(i,j).softCoreLJ.f2(rP2);
					if( rP2 < m_aceCutoff2 ) {
						f_total += aceValue(i,j) * potAce.f2(rP2);
					}
				}
			}
		}

		if( harmonicBasinUse ) {
			for(int i_point = 0; i_point < points2.size(); i_point++) {
				f_total += harmonicBasin.f2(Math::dotSelf(points2(i_point) - points1Center));
			}
		}

		return f_total;

	}


	void gDoubleLoop(const Points& points2,Points& grad2) {

		grad2 = 0;

		for(int i = 0; i < points1.size(); i++) {
			for(int j = 0; j < points2.size(); j++) {
				Point xji = points2(j) - points1(i);
				T_num rP2 = Math::dotSelf(xji);
				if( rP2 < m_cutoffP2 ) {
					grad2(j) += _parNBTableEntry(i,j).softCoreLJ.g2(rP2,xji);
					if( rP2 < m_aceCutoff2 ) {
						grad2(j) += aceValue(i,j) * potAce.g2(rP2,xji);
					}
				}
			}
		}

		if( harmonicBasinUse ) {
			for(int i_point = 0; i_point < points2.size(); i_point++) {
				grad2(i_point) += harmonicBasin.g2(points2(i_point) - points1Center);
			}
		}

	}


protected:

	const ForceParNBTableEntry&
		_parNBTableEntry(int indAtom1,int indAtom2) {
			return fParNBTable(fParAtom1.iType(indAtom1),fParAtom2.iType(indAtom2));
	}


	T_num aceValue(int indAtom1,int indAtom2) {
		return m_aceMatr(fParAtom1.aceType(indAtom1),fParAtom2.aceType(indAtom2));
	}

};


#endif // AT_PRODDL_POTENTIALS_TOTAL_H__
