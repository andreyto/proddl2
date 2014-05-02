//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_DOCKING_H__
#  error DOCKING_MOLSTRUCT.HPP MUST BE INCLUDED FROM WITHIN DOCKING.HPP
#endif

#ifndef PRODDL_DOCKING_MOLSTRUCT_H__
#define PRODDL_DOCKING_MOLSTRUCT_H__



// Cache the intermediate data needed to move the transformation into
// the one that acts on the original coordinates of the ligand:
// Take ligand's transformation found for the reference positions of molecules
// (positions created by 'moveIntoReferencePositions()' method,
// and return transformation relative to the original positions of receptor and ligand.
// In other words, if 'tran' applied to the reference positions of a ligand gives
// a prediction of complex coordinates, then return value applied to the original
// positions (i.e. those passed to the 'init()' method) will give the same complex,
// only rotated/translated as a whole.
// A separate class is created in order to avoid recomputing the inverse of the tR
// on each call.


class ToOriginalFrameTransformer {

public:

	ToOriginalFrameTransformer() {}

	ToOriginalFrameTransformer(const RotationTranslation tR, const RotationTranslation tL):
	m_tR(tR.inverse()),
		m_tL(tL)
	{}

	RotationTranslation operator()(const RotationTranslation& t) const {
		return m_tR * t * m_tL;
	}

protected:

	RotationTranslation m_tR;
	RotationTranslation m_tL;

};


class MolStruct {


public:


	MolStruct(const MolForceParams& params) {

		ATLOG_TRACE_3;

		for(int i_mol = 0; i_mol < N_mol; i_mol++) {	

			trFromIni(i_mol) = RotationTranslation();

		}


		for(int i_mol = 0; i_mol < N_mol; i_mol++) {

			// Positions are an independent copy because we
			// modify them
			pos(i_mol).reference(params.pos(i_mol).copy());

			mass(i_mol).reference(params.mass(i_mol));

		}

		ATLOG_OUT_3("Done copying MolForceParams data");


		minBox = moveIntoReferencePositions();	

	}


	int estimateAngleStepDeg(T_num grid_step) {

		ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

		Geom::Bounding::Diameter<T_num> diamLig(getPosLigand());

		ligSize = diamLig.getSize();

		T_num angDegReal = (grid_step/(ligSize/2))*(180/3.14);

		int angDeg = int(angDegReal);

		if(angDeg < 5)
			angDeg = 5;

		ATLOG_OUT_3(ATLOGVAR(grid_step) << ATLOGVAR(ligSize) << ATLOGVAR(angDeg) << "\n");

		ATLOG_OUT_2("Ideal value of angular step would be: " << angDegReal << \
			". Using instead value " << angDeg << " for performance reasons.");

		return angDeg;

	}



	Points& getPosReceptor() {
		ATLOG_TRACE_3;
		return pos(iRec);

	}

	Points& getPosLigand() {
		ATLOG_TRACE_3;
		return pos(iLig);

	}

	Points& getPos(int i_mol) {
		ATLOG_TRACE_3;
		return pos(i_mol);

	}

	Point centerOfMass(int i_mol) const {
		ATLOG_TRACE_3;

		//return Point(blitz::sum(pos(i_mol) * mass(i_mol))/blitz::sum(mass(i_mol)));
		T_num mass_tot = blitz::sum(mass(i_mol));
		Point x(0.);

		const Points& p = pos(i_mol);
		const Floats& m = mass(i_mol);
		for(int i=pos.lbound(0); i <= pos.ubound(0); ++i) {
			x += p(i)*m(i);
		}
		x /= mass_tot;
		return x;

	}

	const PointPair& getMinBox() const {
		ATLOG_TRACE_3;
		return minBox;

	}

	// Return ligand size (diameter)

	const T_num getSizeLigand() const {
		ATLOG_TRACE_3;
		// probably we forgot to initialize if ligSize == 0

		ATLOG_ASSERT_1(ligSize > 0);

		return ligSize;

	}

	// Take ligand's transformation found for the reference positions of molecules
	// (positions created by 'moveIntoReferencePositions()' method,
	// and return transformation relative to the original positions of receptor and ligand.
	// In other words, if 'tran' applied to the reference positions of a ligand gives
	// a prediction of complex coordinates, then return value applied to the original
	// positions (i.e. those passed to the 'init()' method) will give the same complex,
	// only rotated/translated as a whole.

	RotationTranslation toOriginalFrame(const RotationTranslation& tran) const {
		ATLOG_TRACE_3;
		return  toOrigFrameTransformer(tran);

	}

	// method to get read only access to the underlying object that performs
	// an action of 'toOriginalFrame()'

	const ToOriginalFrameTransformer& getToOriginalFrameTransformer() const {
		ATLOG_TRACE_3;
		return toOrigFrameTransformer;
	}

protected:

	// Move both receptor and ligand into such orientation that the minimum
	// enclosing box is oriented along the coordinate axes, return the lower and
	// upper corners of that box. The cutoff radius of Fft potential 'cutOffFft'
	// is taken into account when calculating the size of the box.
	// Ligand will be additionaly moved to a position where a coordinate origin
	// is at the center of ligand's max diameter.

	PointPair moveIntoMinBox() {

		ATLOG_TRACE_3;

		Geom::Bounding::Box<T_num> boxRec(getPosReceptor());

		RotationTranslation trBox = boxRec.moveToBoundCoordinates();

		for(int i_mol = 0; i_mol < N_mol; i_mol++) {

			trBox(getPos(i_mol));

			trFromIni(i_mol) = trBox*trFromIni(i_mol);

		}

		Geom::Bounding::Diameter<T_num> diamLig(getPosLigand());

		// move ligand such that a center of its diameter is at the coordinate
		// origin

		Translation trToLigDiamCenter = diamLig.getTransformationToBoundCoordinates();

		trToLigDiamCenter(getPos(iLig));

		trFromIni(iLig) = trToLigDiamCenter * trFromIni(iLig);

		ligSize = diamLig.getSize();

		T_num cutOffFft;

		gOptions.get("cutOffFft",cutOffFft);

		T_num gridStep;

		gOptions.get("gridStep",gridStep);

		T_num recPadding = gridStep*2. + ligSize/2. + cutOffFft;

		PointPair diag = boxRec.getDiagonal();

		diag(1) += recPadding;
		diag(0) -= recPadding;

		return diag;

	}


	void randomizeLigandOrientation() {

		ATLOG_TRACE_3;

		T_num startingSpin;

		gOptions.getdefault("startingSpin",startingSpin,T_num(23.3734432));

		Rotation rotSpin(Point(startingSpin,startingSpin*T_num(3.7877),startingSpin*T_num(-4.513432)));

		Point cm = centerOfMass(iLig);

		RotationTranslation rotSpinCm = Translation(cm)*rotSpin*Translation(-1 * cm);

		rotSpinCm(getPosLigand());

		trFromIni(iLig) = rotSpinCm*trFromIni(iLig);

	}


	PointPair moveIntoReferencePositions() {

		ATLOG_TRACE_3;

		randomizeLigandOrientation();

		PointPair ret = moveIntoMinBox();

		toOrigFrameTransformer = ToOriginalFrameTransformer(trFromIni(iRec),trFromIni(iLig));

		return ret;

	}

protected:

	PointsMolSet pos;

	TransformsMolSet trFromIni;

	FloatsMolSet mass;

	PointPair minBox;

	T_num ligSize;

	ToOriginalFrameTransformer toOrigFrameTransformer;

}; // class MolStruct


typedef boost::shared_ptr<MolStruct> PMolStruct;


#endif // PRODDL_DOCKING_MOLSTRUCT_H__

