//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_DOCKING_H__
#  error DOCKING_MOLFORCE.HPP MUST BE INCLUDED FROM WITHIN DOCKING.HPP
#endif

#ifndef PRODDL_DOCKING_MOLFORCE_H__
#define PRODDL_DOCKING_MOLFORCE_H__


class MolForceParams {

public:


	PointsMolSet pos;

	FloatsMolSet mass;

	FloatsMolSet alpha;

	FloatsMolSet eps;

	FloatsMolSet sigma;

	int mix;

};



class MolForce {

public:

	MolForce() {}

	virtual
		~MolForce() {}

	virtual
		void projectMol(int iMol, const Points& positions, VPRawGrid& grids) {

			ATLOG_TRACE_4;

			int iGridFirst = 0;

			if( iMol == iRec ) {
				projectReceptor(positions,grids,iGridFirst);
				//See note on projection with minDist instead of
				//cropping.
				//cropProjectionValues(grids,iGridFirst);
			}
			else {
				projectLigand(positions,grids,iGridFirst);
			}
	}

	// Collect output correlation values from gridsOut(iGridFirst...) 
	// into the total potential values in the gridTot.
	// Implementation must keep in mind that gridTot might be an alias to
	// the one of the gridsOut elements. The latter condition will
	// be always true for the leaf (one-component) MolForce object,
	// so such objects should never copy data into gridTot, although they are allowed
	// to modify 

	void collectTotal(VPRawGrid& gridsOut, PRawGrid gridTot) {

		int iGridFirst = 0;

		collectTotal(gridsOut,iGridFirst,gridTot);

	}


	virtual
		int nGrids() const = 0;


	// The next methods are made public only for implementation
	// purposes (so that MolForceComp class can access them), and
	// should not be used directly by user code. 'projectMol()' must
	// be called instead.

public:


	// projectXXX() must start projecting from grids(iGridFirst) and return
	// the index of the next unused grid (to be used as iGridFirst for the next MolForce object)

	virtual
		int projectReceptor(const Points& positions, VPRawGrid& grids, int iGridFirst) = 0;

	virtual
		int projectLigand(const Points& positions, VPRawGrid& grids, int iGridFirst) = 0;

	virtual
		int collectTotal(VPRawGrid& gridsOut, int iGridFirst, PRawGrid gridTot) = 0;

protected:

	// After we introduced potentials calculated as a sum of several 
	// FFT correlations, it is no longer valid to crop individual grids
	// after projection. For example, cropping in such way would unbalance
	// total LJ potential calculated as a sum of repulsion and attraction grids.
	// Instead, we call projection operator with a minDist parameter equal for
	// both repulsion and attraction grids, thus setting both value to corresponding
	// F(minDist) for all r < minDist.
	// NOTE: when T_num is 'float', FFT is numerically unstable at alpha=0
	// (it produces values strongly dependent on the minDist parameter, both
	// for multicomponent LJ and averaged one-component LJ).
	// Therefore, it is safer to compile this for T_num=double.


	//   virtual 
	//   void cropProjectionValues(VPRawGrid& grids, int iGridFirst) {

	//     ATLOG_TRACE_4;

	//     // TODO: How to find the best cropping limit, where to do the cropping:
	//     // here or somewhere else?

	//     // Crop max repulsion values to avoid overflow ( espesially 
	//     // when T_num is 'float', also happens for large molecules
	//     // even when T_num is 'double' ).

	//     int n_grids = this->nGrids();

	//     for(int i = iGridFirst; i < iGridFirst + n_grids; i++) {

	//       GridArray& arr = grids(i)->getGridArray();

	//       T_num maxLimit = 1.e6;

	//       arr = blitz::where(arr > maxLimit, maxLimit, arr);	

	//       T_num minLimit = - 1.e6;

	//       arr = blitz::where(arr < minLimit, minLimit, arr);

	//     }

	//   }

};


typedef boost::shared_ptr<MolForce> PMolForce;

typedef typename common_types::num_vector_type<PMolForce>::Type VPMolForce;


// Base class for the composition of MolForce objects,
// to implement potentials that can be represented as a sum
// of several components

class MolForceComp : public MolForce {

public:

	// Default implementation - just add components together. In this approach,
	// if different weights are required, they must be used at the projection stage.
	// Child classes can override this method to do some non-additive processing.

	virtual
		int collectTotal(VPRawGrid& gridsOut, int iGridFirst, PRawGrid gridTot) {

			ATLOG_TRACE_4;

			for(int i = 0; i < mfs.size(); i++) {

				iGridFirst = mfs(i)->collectTotal(gridsOut,iGridFirst,gridTot);

			}

			return iGridFirst;

	}

	virtual
		int nGrids() const {

			ATLOG_TRACE_4;

			int n = 0;

			for(int i = 0; i < mfs.size(); i++) {

				n += mfs(i)->nGrids();

			}

			return n;

	}

protected:

	virtual
		int projectReceptor(const Points& positions, VPRawGrid& grids, int iGridFirst) {

			ATLOG_TRACE_4;

			for(int i = 0; i < mfs.size(); i++) {

				iGridFirst = mfs(i)->projectReceptor(positions,grids,iGridFirst);

			}

			return iGridFirst;

	}

	virtual
		int projectLigand(const Points& positions, VPRawGrid& grids, int iGridFirst) {

			ATLOG_TRACE_4;

			for(int i = 0; i < mfs.size(); i++) {

				iGridFirst = mfs(i)->projectLigand(positions,grids,iGridFirst);

			}

			return iGridFirst;

	}

protected:

	VPMolForce mfs;

}; // class MolForceMulti


// Leaf class of MolForce tree

class MolForceSingle : public MolForce {

public:

	virtual
		int collectTotal(VPRawGrid& gridsOut, int iGridFirst, PRawGrid gridTot) {

			ATLOG_TRACE_4;

			PRawGrid gridOut = gridsOut(iGridFirst);

			ATLOG_OUT_4(ATLOGVAR(iGridFirst) << ATLOGVAR(gridOut) << ATLOGVAR(gridTot));

			// Never add to itself

			if( gridOut != gridTot ) {

				// That might also add the nonsense FFT padding together, but it should not do any harm

				gridTot->getGridArray() += gridOut->getGridArray();

			}

			return iGridFirst + 1;

	}


	virtual
		int nGrids() const {

			ATLOG_TRACE_4;

			return 1;

	}

};


// Implements Receptor projection of radial field.
// Ligand projection still has to be implemented by children of this class

template<class TRadialFieldR>
class MolForceRadialR : public MolForceSingle {

public:

	MolForceRadialR():
	  minDist(-1)
	  {}

protected:

	virtual
		int projectReceptor(const Points& positions, VPRawGrid& grids, int iGridFirst) {

			ATLOG_TRACE_4;

			Grid& grid = *grids(iGridFirst);

			//TODO: introduce external parameter for minDist

			T_num _minDist = minDist;
			if( _minDist < 0 ) {
				_minDist = grid.minSpatialStep() / 4.;
			}

			ATLOG_OUT_4(ATLOGVAR(_minDist));

			projectorR.init(grid,_minDist);

			projectorR.projectFields(fieldsR,positions);

			return iGridFirst + 1;

	}


protected:

	typedef TRadialFieldR RadialFieldR;

	typedef typename common_types::num_vector_type<RadialFieldR>::Type RadialFieldsR;

	RadialFieldsR fieldsR;

	typedef typename Projector::template ProjectorRadialField<RadialFieldR> ProjectorR;

	ProjectorR projectorR;

	T_num minDist;

}; // class MolForceRadialR


// Implements Ligand projection as addition of some real number value
// defined for each atom to the nearest grid point

template<class TRadialFieldR>
class MolForceRadial : public MolForceRadialR<TRadialFieldR> {

public:


protected:

	virtual
		int projectLigand(const Points& positions, VPRawGrid& grids, int iGridFirst) {

			ATLOG_TRACE_4;

			Grid& grid = *grids(iGridFirst);

			grid = 0.;

			for(int i_atom = 0; i_atom < positions.size(); i_atom++) {

				grid(positions(i_atom)) += fieldsL(i_atom);

			}

			return iGridFirst + 1;

	}

protected:

	Floats fieldsL;

}; // class MolForceRadial


class MolForceLJ : public MolForceRadial<LjPot> {

public:

    typedef MolForceRadial<LjPot> Base; 

	MolForceLJ(const MolForceParams& params) {

		ATLOG_TRACE_4;

		int mix = params.mix;

		const Floats& eps_r = params.eps(iRec);
		const Floats& sigma_r = params.sigma(iRec);
		const Floats& alpha_r = params.alpha(iRec);
		// find average values of LJ params for ligand molecule
		// TODO: possibly, surface only atoms should be used here

		T_num eps_avg_l = blitz::mean(params.eps(iLig));
		T_num sigma_avg_l = blitz::mean(params.sigma(iLig));

		T_num cutOffFft;

		gOptions.get("cutOffFft",cutOffFft);	  


        Base::fieldsR.resize(eps_r.size());


		for(int i_atom = 0; i_atom < eps_r.size(); i_atom++) {

			T_num sigma;

			if( mix == PotentialsT::LJ_MIX_0 ) {
				sigma = (sigma_r(i_atom)+sigma_avg_l)/2.;
			}
			else {
				sigma = std::sqrt(sigma_r(i_atom)*sigma_avg_l);
			}

            Base::fieldsR(i_atom) = LjPot(
				sigma, 
				std::sqrt(eps_r(i_atom)*eps_avg_l), 
				alpha_r(i_atom),
				cutOffFft
				);

			////Made results slightly worse
			// 	  fieldsR(i_atom) = LjPot(
			// 				    (sigma_r(i_atom)+sigma_avg_l)/2., 
			// 				    std::sqrt(eps_r(i_atom)), 
			// 				    alpha_r(i_atom),
			// 				    cutOffFft
			// 				    );

		}


		const Floats& eps_l = params.eps(iLig);

        Base::fieldsL.resize(eps_l.size());

		for(int i_atom = 0; i_atom < eps_l.size(); i_atom++) {

			////Made results sligthly worse, see above for receptor
			// 	  fieldsL(i_atom) = std::sqrt(eps_l(i_atom));

            Base::fieldsL(i_atom) = 1;

		}

	}

}; // class MolForceLJ


// Repulsive part of LJ

class MolForceLJRep : public MolForceRadial<LjPotRep> {

public:

    typedef MolForceRadial<LjPotRep> Base;

	MolForceLJRep(const MolForceParams& params) {

		ATLOG_TRACE_4;

		int mix = params.mix;

		ATLOG_ASSERT_1(mix == PotentialsT::LJ_MIX_1);

		const Floats& eps_r = params.eps(iRec);
		const Floats& sigma_r = params.sigma(iRec);
		const Floats& alpha_r = params.alpha(iRec);
		// find average values of LJ params for ligand molecule
		// TODO: possibly, surface only atoms should be used here

		T_num sigma_avg_l = blitz::mean(params.sigma(iLig));

		T_num cutOffFft;

		gOptions.get("cutOffFft",cutOffFft);	  


        Base::fieldsR.resize(eps_r.size());


		for(int i_atom = 0; i_atom < eps_r.size(); i_atom++) {

			T_num sigmaAlpha = std::sqrt(sigma_r(i_atom)*sigma_avg_l);

            Base::fieldsR(i_atom) = LjPotRep(
				std::sqrt(sigma_r(i_atom)),
				sigmaAlpha,
				std::sqrt(eps_r(i_atom)), 
				alpha_r(i_atom),
				cutOffFft
				);

		}


		const Floats& eps_l = params.eps(iLig);
		const Floats& sigma_l = params.sigma(iLig);

        Base::fieldsL.resize(eps_l.size());

		for(int i_atom = 0; i_atom < eps_l.size(); i_atom++) {

			// pow6(sigma) == pow12(sqrt(sigma))

            Base::fieldsL(i_atom) = std::sqrt(eps_l(i_atom)) * Math::pow6(sigma_l(i_atom));

		}

	}

}; // class MolForceLJRep



// Attractive part of LJ

class MolForceLJAttr : public MolForceRadial<LjPotAttr> {

public:

    typedef MolForceRadial<LjPotAttr> Base;

	MolForceLJAttr(const MolForceParams& params) {

		ATLOG_TRACE_4;

		int mix = params.mix;

		ATLOG_ASSERT_1(mix == PotentialsT::LJ_MIX_1);

		const Floats& eps_r = params.eps(iRec);
		const Floats& sigma_r = params.sigma(iRec);
		const Floats& alpha_r = params.alpha(iRec);
		// find average values of LJ params for ligand molecule
		// TODO: possibly, surface only atoms should be used here

		T_num sigma_avg_l = blitz::mean(params.sigma(iLig));

		T_num cutOffFft;

		gOptions.get("cutOffFft",cutOffFft);	  


        Base::fieldsR.resize(eps_r.size());


		for(int i_atom = 0; i_atom < eps_r.size(); i_atom++) {

			T_num sigmaAlpha = std::sqrt(sigma_r(i_atom)*sigma_avg_l);

            Base::fieldsR(i_atom) = LjPotAttr(
				std::sqrt(sigma_r(i_atom)),
				sigmaAlpha,
				std::sqrt(eps_r(i_atom)), 
				alpha_r(i_atom),
				cutOffFft
				);

		}


		const Floats& eps_l = params.eps(iLig);
		const Floats& sigma_l = params.sigma(iLig);

        Base::fieldsL.resize(eps_l.size());

		for(int i_atom = 0; i_atom < eps_l.size(); i_atom++) {

			// pow3(sigma) == pow6(sqrt(sigma))

            Base::fieldsL(i_atom) = std::sqrt(eps_l(i_atom)) * Math::pow3(sigma_l(i_atom));

		}

	}

}; // class MolForceLJAttr




class MolForceLJComp : public MolForceComp {

public:

	MolForceLJComp(const MolForceParams& params) {

		ATLOG_TRACE_4;

		this->mfs.resize(2);

		this->mfs(0).reset(new MolForceLJRep(params));
		this->mfs(1).reset(new MolForceLJAttr(params));

	}


}; // class MolForceLJComp


// Potential to test correctness of implementation:
// Just adds two copies of a full LJ potential -
// Must provide same coordinates in results file
// but with doubled energy values.

class MolForceLJCompTest : public MolForceComp {

public:

	MolForceLJCompTest(const MolForceParams& params) {

		ATLOG_TRACE_4;

		this->mfs.resize(2);

		this->mfs(0).reset(new MolForceLJ(params));
		this->mfs(1).reset(new MolForceLJ(params));

	}


}; // class MolForceLJComp




#endif // PRODDL_DOCKING_MOLFORCE_H__

