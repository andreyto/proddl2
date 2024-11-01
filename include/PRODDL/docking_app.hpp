//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_DOCKING_H__
#  error DOCKING_APP.HPP MUST BE INCLUDED FROM WITHIN DOCKING.HPP
#endif

#ifndef PRODDL_DOCKING_APP_H__
#define PRODDL_DOCKING_APP_H__

class App {

public:

public:

	App() {

	}

	virtual ~App() {

		ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

	}

public:

	void init(const MolForceParams& mfParams) {

		ATLOG_TRACE_3;

		int runTimeLogLevel;

		gOptions.getdefault("logLevel",runTimeLogLevel,ATLOG_LEVEL_1);

		Logger::setRunTimeLevel(runTimeLogLevel);

		gOptions.getdefault("testMode",testMode,0);

		ATLOG_OUT_1(ATLOGVAR(ATLOG_LEVEL) << ATLOGVAR(Logger::getRunTimeLevel()) << ATLOGVAR(testMode));

		pmolStruct.reset(new MolStruct(mfParams));

		ATLOG_OUT_4("pmolStruct initialized"); 

		std::string potentialName;

		gOptions.getdefault("potentialName",potentialName,"ljComp");

		MolForce * p_mf = 0;

		if( potentialName == "ljComp" ) {
			p_mf = new MolForceLJComp(mfParams);
		}
		else if( potentialName == "ljAvg" ) {
			p_mf = new MolForceLJ(mfParams);
		}
		else {
			throw not_supported_error("Unknown 'potentialName' parameter value: " + potentialName);
		}

		pmolForce.reset(p_mf);

		ATLOG_OUT_2("Using potential defined on " << pmolForce->nGrids() << " grids.");

		gOptions.getdefault("gridStep",gridStep,T_num(1.5));

		gOptions.getdefault("angleStepDeg",angleStepDeg,-1);

		if( angleStepDeg < 0 ) {

			angleStepDeg = pmolStruct->estimateAngleStepDeg(gridStep);

		}

		ATLOG_OUT_3(ATLOGVAR(angleStepDeg));

		gOptions.getdefault("maxNTrans",maxNTrans,10);

	}

	virtual void run() = 0;

	virtual bool isForeman() const = 0;


protected:

	PMolStruct pmolStruct;

	PMolForce pmolForce;

	T_num gridStep;

	//angle step in degrees, rounded to whole degrees
	int angleStepDeg;

protected:

	// Maximum number of translations to select in worker and send back to the foreman process.

	int maxNTrans;

	// If not 0, this switch puts the programm into test mode

	int testMode;

};


class Foreman : public App {

public:

	typedef App Base;
	typedef Foreman Self;

	typedef RotFftScanIO_Collector<RotFftScanIO_Bin> RotFftScanIO_CollectorT;

protected:


	//TODO: std::binary_function: pass args by value?
	struct cmp_rt_val : public std::binary_function<RotTranValue,RotTranValue,bool> {

		cmp_rt_val()
		{}

		bool operator()(const RotTranValue& x,const RotTranValue& y) {

			return (x.value < y.value);

		}

		bool operator()(const RotTranValue& x,T_num y) {

			return (x.value < y);

		}

	}; // struct cmp_rt_val

	typedef BoundPriorityQueue<RotTranValue,std::vector<RotTranValue>,cmp_rt_val,T_num,cmp_rt_val> RotTranValQueueType;

public:

	typedef typename RotTranValQueueType::container_type ResultsContainerType;

public:

	void init(const MolForceParams& mfParams) {

		ATLOG_TRACE_3;

		Base::init(mfParams);


		tranValues.resize(this->maxNTrans);

		gOptions.getdefault("maxNRigidMatches",maxNRigidMatches,20000);

		rtvalQueue.init(maxNRigidMatches,T_num(1e16));

		resultsInOriginalFrame = false;

	}

	void run() {

		ATLOG_TRACE_3;    

		ATLOG_STD_EXCEPTIONS_TRY();

		startCollectRotTranVals();

		while(nextCollectRotTranVals()) {


			int nTranValues = tranValues.size();

			ATLOG_OUT_4("Received results for one rotation: " <<  \
				ATLOGVAR(nTranValues));

			for(int iTranVal = 0; iTranVal < nTranValues; iTranVal++) {

				const TranValue& tranVal = tranValues(iTranVal);

				RotTranValue rtVal;

				rtVal.tran = tranVal.tran * doneRot;

				rtVal.value = tranVal.value;

				rtvalQueue.push(rtVal);

				// DEBUG:
				if( iTranVal >= ( nTranValues - 4 ) ) {

					ATLOG_OUT_4(ATLOGVAR(iTranVal) \
						<< ATLOGVAR(rtVal.tran) \
						<< ATLOGVAR(rtVal.value));

				}

			}
		}

		finishCollectRotTranVals();

		rtvalQueue.sort();

		ATLOG_STD_EXCEPTIONS_CATCH();

	}

	bool isForeman() const {

		return true;

	}

	void startCollectRotTranVals() {
		ATLOG_TRACE_3;
		std::string fft_rot_scan_list;
		gOptions.get("fft_rot_scan_list",fft_rot_scan_list);
		ATALWAYS(m_io_rot_scan_coll.get() == 0,"Input object already exists");
		m_io_rot_scan_coll.reset(new RotFftScanIO_CollectorT(fft_rot_scan_list));
	}

	bool nextCollectRotTranVals() {
		ATLOG_TRACE_3;
		ATALWAYS(m_io_rot_scan_coll.get(),"Input object does not exist");
		return m_io_rot_scan_coll->read_record(doneRot,tranValues);
	}

	void finishCollectRotTranVals() {
		ATLOG_TRACE_3;
		ATALWAYS(m_io_rot_scan_coll.get(),"Input object does not exist");
		m_io_rot_scan_coll.reset();
	}

	const ResultsContainerType& getResults() {
		ATLOG_TRACE_3;

		ResultsContainerType& results = rtvalQueue.data();

		if ( ! resultsInOriginalFrame ) {

			for(int i = 0; i < results.size(); i++) {

				results[i].tran = this->pmolStruct->toOriginalFrame(results[i].tran);

			}

			resultsInOriginalFrame = true;

		}

		return results;

	}

	void
		writeResults(const std::string& fileName, char format) {
			ATLOG_TRACE_3;

			const ResultsContainerType& results = getResults();

			IORigid<T_num> io;

			io.writeCoords(fileName,results.begin(),results.size(),format);
			//io.readCoords(fileName,format);
			//io.writeCoords(fileName,format);

	}


protected:


	// restricted size priority queue to accumulate transformations with corresponding
	// values of a target function

	RotTranValQueueType rtvalQueue;

	// flag is true when transformations in rtvalQueue were already converted
	// into original coordinate frame (that should be done only once)

	bool resultsInOriginalFrame;

	// buffer to receive processed ligand rotation 

	Rotation doneRot;

	// buffer to receive translations found for a given ligand rotation

	TranValues tranValues;

	// maximum number of matches to select during rigid body docking

	int maxNRigidMatches;

	/// Reader for a collection of rotation scan result files

	boost::scoped_ptr<RotFftScanIO_CollectorT> m_io_rot_scan_coll;

};


class Worker : public App {

public:

	typedef App Base;
	typedef Worker Self;


public:

	void init(const MolForceParams& mfParams) {

		ATLOG_TRACE_3;

		Base::init(mfParams);

		std::string anglesFile;
		gOptions.get("anglesFile",anglesFile);

		rotGrid.reference(Geom::loadRotationalGrid<T_num>(anglesFile));

		rotToRun.clear();

		int fft_rot_grid_start = 0;
		gOptions.get("fft_rot_grid_start",fft_rot_grid_start);

		ATALWAYS(fft_rot_grid_start >= 0 && fft_rot_grid_start < rotGrid.size(),\
			"Rotation start index is out of bound");

		int fft_rot_grid_end = 0;
		gOptions.get("fft_rot_grid_end",fft_rot_grid_end);

		ATALWAYS(fft_rot_grid_end >= fft_rot_grid_start,\
			"rotation end index is out of bound");

		 if(fft_rot_grid_end > rotGrid.size()) {
			 fft_rot_grid_end = rotGrid.size();
		 }

		for( int i = fft_rot_grid_start; i < fft_rot_grid_end; i++) 
			rotToRun.push_back(rotGrid(i));


		//TODO: some intelligent estimate for default
		//values of maxNTrans and maxValCorr

		T_num maxValCorr;
		gOptions.getdefault("maxValCorr",maxValCorr,T_num(0));
		pfft.reset(new FFTCorrelators(this->pmolStruct->getMinBox(),this->gridStep,this->pmolForce->nGrids()));

		prepareReceptor();
		Grid& corrGrid = *(pfft->getGridTot());

		bool doClusterTranslations;
		gOptions.getdefault("doClusterTranslations",doClusterTranslations,false);

		int maxNTransInp;

		gOptions.getdefault("maxNTransInp",maxNTransInp,this->maxNTrans*100);

		if( doClusterTranslations ) {

			T_num clusterRadiusTrans;
			gOptions.getdefault("clusterRadiusTrans",clusterRadiusTrans,5.0);

			pTranClust.reset(new TranClust(maxNTransInp, this->maxNTrans, clusterRadiusTrans));

		}

		int nMultimer;
		gOptions.getdefault("nMultimer",nMultimer,0);

		if( nMultimer > 1 ) {

			T_num maxRmsdSymm;

			gOptions.getdefault("maxRmsdSymm",maxRmsdSymm,8.0);

			pTranSymm.reset(new TranSymm(nMultimer, maxRmsdSymm, this->pmolStruct->getSizeLigand()/2,
				this->pmolStruct->getToOriginalFrameTransformer()));

		}


		fftProc.init(corrGrid, maxNTransInp, maxValCorr);

		resetLigPosRot();

	}

	void prepareReceptor() {

		ATLOG_TRACE_3;

		VPRawGrid& recGrid = pfft->getGridsRec();
		this->pmolForce->projectMol(iRec,this->pmolStruct->getPosReceptor(),recGrid);

		pfft->preprocessReceptor();

	}

	void findTranslationalMinima() {

		ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

		VPRawGrid& ligGrid = pfft->getGridsLig();

		this->pmolForce->projectMol(iLig,ligPosRot,ligGrid);

		pfft->correlate();

		this->pmolForce->collectTotal(pfft->getGridsOut(),pfft->getGridTot());

		fftProc.selectFromFFT();

		if( pTranSymm ) {

			pTranSymm->init(currentRot);

			fftProc.postProcess(pTranSymm);

		}

		if( pTranClust ) {

			fftProc.postProcess(pTranClust);

		}

		//dbg::out(dbg::info) << dbg::indent() << ATLOGVAR(ligGrid.getGridArray()) << '\n';

	}


	void scanRotation() {

		ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

		resetLigPosRot();

		currentRot(ligPosRot);

		findTranslationalMinima();

		TranValues tranVals = fftProc.getTranValuesFilled(this->maxNTrans);

		ATLOG_OUT_4("Output translations: " << ATLOGVAR(tranVals.size()));

		outputScannedRotation(currentRot,tranVals);

	}

	void outputScannedRotation(const Rotation& rot, const TranValues& tranVals) {
		ATLOG_TRACE_3;
		ATALWAYS(m_io_rot_scan->write_record(rot,tranVals),"Output failed");
	}

	void startOutputScannedRotations() {
		ATLOG_TRACE_3;
		std::string file_name;
		gOptions.get("fft_rot_scan_res",file_name);
		ATALWAYS(m_io_rot_scan.get() == 0,"Output object already exists");
		m_io_rot_scan.reset(new RotFftScanIO_Bin(file_name,std::fstream::out));
	}

	void finishOutputScannedRotations() {
		ATLOG_TRACE_3;
		ATALWAYS(m_io_rot_scan.get(),"Output object does not exist");
		m_io_rot_scan.reset();
	}

	void run() {

		ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

		ATLOG_STD_EXCEPTIONS_TRY();

		startOutputScannedRotations();

		//process all rotations
		while(! rotToRun.empty()) {

			currentRot = rotToRun.front();

			ATLOG_OUT_3(ATLOGVAR(currentRot.eulerAngles()));
			scanRotation();

			rotToRun.pop_front();

			ATLOG_OUT_3(ATLOGVAR(rotToRun.size()));

		}

		finishOutputScannedRotations();

		ATLOG_STD_EXCEPTIONS_CATCH();

	}


	bool isForeman() const {

		return false;

	}


	// This method must be used for testing only if
	// called from Python, because I did not have time to figure out
	// the life time and copy semantics of returned 
	// smart pointer in Boost Python library.

	PFFTCorrelator getFFTCorrelator(int iFft) {
		return pfft->getFfts()(iFft);
	}

	// This testing function will leave projections
	// of receptor and ligand in corresponding grids
	// of 'fft' object.

	void testProjection() {

		ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

		VPRawGrid& recGrid = pfft->getGridsRec();
		this->pmolForce->projectMol(iRec,this->pmolStruct->getPosReceptor(),recGrid);

		VPRawGrid& ligGrid = pfft->getGridsLig();
		this->pmolForce->projectMol(iLig,ligPosRot,ligGrid);

		//DEBUG:
		//Checking for overflow. If everything is fine,
		//correlation with a single 1 point should  give
		//the orginal receptor grid.
		//GridArray& recArr = recGrid.getGridArray();
		//recArr = blitz::where(recArr > 100,100,recArr);
		//ligGrid = 0;
		//ligGrid(Point(0,0,0)) = 1;

	}


protected:

	void resetLigPosRot() {
		ATLOG_TRACE_3;
		ligPosRot.reference(this->pmolStruct->getPosLigand().copy());
	}

protected:

	// all rotations, forming a grid with a defined angle step
	Rotations rotGrid;


	// rotations to process during docking

	DequeRotations rotToRun;

	PFFTCorrelators pfft;

	CorrelationProcessor fftProc;

	PTranProcessor pTranClust;

	PTranSymm pTranSymm;

	// ligand atomic coordinates for the current rotation

	Points ligPosRot;


	// current rotation as Rotation object

	Rotation currentRot;

	/// IO object for IPC

	boost::scoped_ptr<RotFftScanIO_Bin> m_io_rot_scan;


};

#endif // PRODDL_DOCKING_APP_H__

