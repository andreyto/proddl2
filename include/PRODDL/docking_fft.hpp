//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_DOCKING_H__
#  error DOCKING_FFT.HPP MUST BE INCLUDED FROM WITHIN DOCKING.HPP
#endif

#ifndef PRODDL_DOCKING_FFT_H__
#define PRODDL_DOCKING_FFT_H__


class TranProcessor {

public:

  virtual
  void process(int nInp, TranValues& tranVals, int& nOut) = 0;
  virtual ~TranProcessor() {}
}; // class TranProcessor


typedef boost::shared_ptr<TranProcessor> PTranProcessor;


class TranClust : public TranProcessor, public boost::noncopyable {

protected:

  typedef ClusterMatrix<T_num,N_dim> ClusterTrans;

  typedef typename ClusterTrans::Matrix Matrix;


public:


  TranClust(int maxInpN, int maxOutN, T_num clusterRadius) {

    matrTran.resize(maxInpN,N_dim);

    weightsTran.resize(maxInpN);

    indRowTran.resize(maxInpN);

    clustDens.resize(maxInpN);

    clustRad = clusterRadius;

    nOutMax = maxOutN;

  }

  virtual
  void process(int nInp, TranValues& tranVals, int& nOut) {

    ATLOG_ASSERT_1(nInp <= tranVals.rows() && nInp <= matrTran.rows());

      // view matrix to be clustered as a vector of points
      Points pointsTran(blitz::viewWithFoldedComponent<Point>(matrTran));

      for(int i = 0; i < nInp; i++) {

	const TranValue& tranVal = tranVals(i);

	pointsTran(i) = tranVal.tran.getVector();
	weightsTran(i) = - tranVal.value;

      }

      // Create smaller views into clustering data

      Matrix m_cl(matrTran(blitz::Range(0,nInp-1),blitz::Range::all()));
      Floats w_cl(weightsTran(blitz::Range(0,nInp-1)));

      clustTran.cluster(m_cl,w_cl,clustRad);

      clustTran.getClusters(indRowTran,clustDens);

      // we rely on the property that results of getClusters
      // are already sorted

      nOut = clustTran.numClusters();

      ATLOG_OUT_3("Found " << nOut << " clusters among " << nInp << " translations.");
      
      if( nOut > tranVals.rows() ) {
	nOut = tranVals.rows();
      }

      if( nOut > nOutMax ) {
	nOut = nOutMax;
      }

      for( int i = 0; i < nOut; i++ ) {
	
	TranValue& tranVal = tranVals(i);
	int indRow = indRowTran(i);
	tranVal.tran = Translation(pointsTran(indRow));
	tranVal.value = - weightsTran(indRow);

      }

  }

protected:

  // radius to use for clustering

  T_num clustRad;

  // clustering object

  ClusterTrans clustTran;

  // intermediate input matrix N_queuexN_dim for clustering call

  Matrix matrTran;

  // intermediate input array for weight for clustering call

  Floats weightsTran;

  // intermediate output array for indexes of cluster leaders

  Ints indRowTran;

  // intermediate output array for cluster density

  Floats clustDens;

  // maximum number of output elements

  int nOutMax;

}; // class TranClust



// Filter correlation results based on C(n) symmetry constraint

class TranSymm : public TranProcessor, public boost::noncopyable {

protected:


public:


  TranSymm(int nSymm, T_num maxRmsdSymm, T_num molRadius, 
	   ToOriginalFrameTransformer toOriginalFrameTransformer):
    m_maxRmsdSymm(maxRmsdSymm),
    m_symmChecker(nSymm,molRadius),
    m_toOriginalFrameTransformer(toOriginalFrameTransformer)
  {
    // m_rot will be set to an Identity
    // from its default ctor
  }

  
  void
  init(Rotation rot) {
    m_rot = rot;
  }


  virtual
  void process(int nInp, TranValues& tranVals, int& nOut) {

    ATLOG_ASSERT_1(nInp <= tranVals.rows());

    ATLOG_OUT_3("Starting " << m_symmChecker.getNSymm() << "-fold symmetry filtering." << ATLOGVAR(m_maxRmsdSymm));
      
    int i_out = 0;

    for( int i_inp = 0; i_inp < nInp; i_inp++ ) {
	
      const TranValue& tranValInp = tranVals(i_inp);

      // translation applied after rotation
      RotationTranslation transform = m_toOriginalFrameTransformer(tranValInp.tran * m_rot);

      T_num rmsd = m_symmChecker(transform);

      if( rmsd < m_maxRmsdSymm ) {

	TranValue& tranValOut = tranVals(i_out);
	tranValOut = tranValInp;
	i_out++;

      }

    }

    nOut = i_out;

    ATLOG_OUT_2("Post-FFT symmetry filter retained " << nOut << " translations out of " << nInp);

  }

protected:

  // Rotation used for this FFT

  Rotation m_rot;


  // Maximum allowed rmsd from symmetry

  T_num m_maxRmsdSymm;

  Geom::SymmetryCheckerCn<T_num> m_symmChecker;

  // Copy of an object created by MolStruct, that can be
  // used to convert each transformation into the one
  // that is relative to the initial positions of
  // the receptor and ligand, which should be identical
  // for the homo-multimer (the positions that are used
  // during FFT are not, because the receptor is moved
  // into the minimum enclosing box frame, and the ligand's
  // orientation is randomized). The symmetry test used here is
  // only applicable to the original equal positions.

  ToOriginalFrameTransformer m_toOriginalFrameTransformer;

}; // class TranSymm


typedef boost::shared_ptr<TranSymm> PTranSymm;


class CorrelationProcessor {

protected:

  struct cmp : public std::binary_function<int,int,bool> {

    cmp() : 
      data(0) 
    {}

    cmp(const T_num *_data) :
      data(_data)
    {}

    bool operator()(int x,int y) {

      return (data[x] < data[y]);

    }

    // we will use this class as a _LimitCompare template argument to the BoundPriorityQueue also
    bool operator()(int x,T_num y) {

      return (data[x] < y);

    }

  protected:

    const T_num *data;

  }; // struct cmp

  typedef BoundPriorityQueue<int,std::vector<int>,cmp,T_num,cmp> QueueType;

  typedef Math::WrappedIndex<T_num,N_dim> WrappedIndexType;


public:


  CorrelationProcessor()
  {
    ATLOG_TRACE_3;
    clear();
  }

  void init(Grid& grid, int maxInpN, T_num maxInpVal) {

    ATLOG_TRACE_3;

    clear();

    p_grid = &grid;

    gridRawData = grid.getGridArray().dataFirst();

    ATALWAYS(blitz::all(grid.getGridArray().lbound() == 0),"Only zero-bazed arrays are supported");
    //TODO: either make a check that the array is in C storage order,
    // or generalize corresponding array functions (such as index computation).

    ATLOG_OUT_3(ATLOGVAR(maxInpN) << ATLOGVAR(maxInpVal));

    queue.init(maxInpN,maxInpVal,cmp(gridRawData),cmp(gridRawData));

    //Important: shape of the grid, as returned by grid.getLogicalShape()
    //can be smaller (and will be, in case of in-place transforms) than the shape
    //of the grid's internal array (because of the required padding, for instance).
    //We assume here that the grid's shape is the size of FFT transform.

    wrappedIndex = WrappedIndexType(grid.getLogicalShape());

    stride = grid.getGridArray().stride();

    corrTranResults.resize(maxInpN);

  }



  void
  selectFromFFT(bool do_sort=true) {

    ATLOG_TRACE_4;

    selectIntoQueue();

    if( do_sort )
      queue.sort();

    ATALWAYS(corrTranResults.isStorageContiguous(),"");

    int n_queue = queue.data().size();


    for(int i = 0; i < n_queue; i++) {

      int ind = queue.data()[i];

      TranValue& corrTranResult = corrTranResults(i);

      corrTranResult.tran = Translation(offsetToSpaceCoords(ind,corrTranResult.value));

    }

    nOut = n_queue;

    ATLOG_OUT_4(ATLOGVAR(nOut) << ATLOGVAR(queue.getSizeLimit()));

  }

  // Apply the 'pproc' object to post-prosess correlation values
  // stored in 'corrTranResults'. After the call,
  // the final results are stored at the beginning of 'corrTranResults',
  // (there are might be less than before the call),
  // and the this->size() == nOut output argument to pproc->process(...).
  // Therefore, this method can be called again with another 'pproc' object,
  // thus applying a chain of postprocessing objects.

  void postProcess(PTranProcessor pproc) {
    ATLOG_TRACE_4;

    pproc->process(this->size(),corrTranResults,nOut);

  }

  int rawLogicalCoordsToOffset(const IntPoint& ind) {

    return blitz::dot(stride,ind);

  }

  Point offsetToSpaceCoords(int ind, T_num& val) {

    val = gridRawData[ind];

    IntPoint indN;
    for(int dim = 0; dim < N_dim; dim++) {

      indN(dim) = ind/stride(dim);

      ind %= stride(dim);

    }

    indN = wrappedIndex.unwrap(indN);

    // Logical index represents the displacement of a ligand
    // regardless of the spatial coordinate origin,
    // hence using toSpatialDiff() here.

    return p_grid->getGeometry().toSpatialDiff(indN);

  }

  int size() const {

    ATLOG_TRACE_4;
	
    return nOut;
	
  }

  void clear() {

    ATLOG_TRACE_4;
	
    queue.clear();
    nOut = 0;

  }

  QueueType&
  getQueue() {

    return queue;

  }

  // not the entire returned array may be filled with data:
  // the number of meaningfull data records is 'this->size()'

  TranValues&
  getTranValues() {
    ATLOG_TRACE_4;

    return corrTranResults;
	
  }

  // return view of the output array that contains only records filled with data
  // but no more than nMax records

  TranValues
  getTranValuesFilled(int nMax) {
    ATLOG_TRACE_4;

    int nSize = size();

    if( nSize > nMax) {
      nSize = nMax;
    }

    if (nSize > 0) {
      return TranValues(corrTranResults,blitz::Range(0,nSize-1));
    }
    else {
      TranValues empty;
      return empty;
    }
	
  }

protected:

  void
  selectIntoQueue() {

    ATLOG_TRACE_4;
	
    clear();

    //This code will exclude padded part of the array.

    for(IntPoint ind = wrappedIndex.before_first(); wrappedIndex.next(ind); ) {
	  
      queue.push(rawLogicalCoordsToOffset(ind));
	  
    }

    ATLOG_OUT_4(ATLOGVAR(queue.size()));

  }


protected:

  QueueType queue;

  const Grid *p_grid;

  const T_num *gridRawData;

  WrappedIndexType wrappedIndex;

  IntPoint stride;

  // output array for finally selected translations with values

  TranValues corrTranResults;

  // actual number of finally selected rotations in corrTranResults,
  // shoudl be less or equal than corrTranResults.size()

  int nOut;

};



// The sequence of calls to methods of this class is:
// FFTCorrelator corr;
// corr.init(...)
// Grid& recGrid = corr.getGrid(iGridRec);
// ... project receptor ...
// corr.preprocessReceptor(); // does the direct FFT of receptor
// Grid& ligGrid = corr.getGrid(iGridLig);
// while(rotations) {
//    ... project ligand ...
//    corr.correlate(); // correlation function is in ligGrid
//    ... select best correlation points ...
// }

class FFTCorrelator {

public:

  enum GridIndex { iGridRec, iGridLig };

  enum { N_grids = 2 };

  // ligand grid will hold data for reciprocal transform

  enum { iGridC2R = iGridLig, iGridOut = iGridC2R };

  // small array of Grids

  typedef typename common_types::point_type<Grid,N_grids>::Type GridsSet;

  typedef typename common_types::point_type<GridArrayC,N_grids>::Type GridArraysCSet;

  typedef FftwPlan<T_num> FftwPlanType;

  typedef typename FftwPlanType::T_complex FftwComplex;

  typedef typename common_types::point_type<FftwPlanType,N_grids>::Type FftwPlansSet;

public:

  void init(const PointPair& boxDiag, T_num gridStep) {

    ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

    // find the smallest size for FFT grid

    Point grStep(gridStep);
    typename Grid::Geom gridGeom(grStep,boxDiag(0));


    Math::FFTW_Size::findBestSize<N_dim>(gridGeom.toLogical(boxDiag(1)),
					 fftSize);

    ATLOG_SWITCH_3(dbg::out(dbg::info) << dbg::indent() \
		   << ATLOGVAR(boxDiag) \
		   << ATLOGVAR(grStep) \
		   << ATLOGVAR(fftSize) << "\n");

    gridGeom = typename Grid::Geom(Point(gridStep),fftSize);

    // create the grids

    for( int i_grid = 0; i_grid < N_grids; i_grid++) {

      IntPoint fftComplPhysSize = fftSize;
      fftComplPhysSize(2) = fftSize(2)/2 + 1;
      arraysC(i_grid).reference(GridArrayC(fftComplPhysSize));
      arraysC(i_grid) = 0; // std::complex<T_num>(0,0);
      IntPoint fftRealPhysSize = fftComplPhysSize;
      // see FFTW 3 docs about padding the last dimension of real arrays for
      // in-place real to complex transforms
      fftRealPhysSize(2) *= 2;

      grids(i_grid).init(gridGeom,fftSize,
			 GridArray(reinterpret_cast<T_num*>(arraysC(i_grid).dataFirst()),
				   fftRealPhysSize,blitz::neverDeleteData));

      T_num * pReal = grids(i_grid).getGridArray().dataFirst();
      FftwComplex *pCompl = reinterpret_cast<FftwComplex*>(arraysC(i_grid).dataFirst());

      // grids(i_grid).getGridArray()(IntPoint(0)) = 1;

      ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() << ATLOGVAR(i_grid) \
		     << ATLOGVAR(pReal) << ATLOGVAR(pCompl) \
		     << ATLOGVAR(*pReal) \
		     << ATLOGVAR((*pCompl)[0]) << ATLOGVAR((*pCompl)[1]) << "\n");

      fftwPlansR2C(i_grid).dft_r2c(fftSize.length(),fftSize.dataFirst(),
				   pReal,
				   pCompl,
				   FFTW_PATIENT);

    }


    T_num * pReal = grids(iGridC2R).getGridArray().dataFirst();
    FftwComplex *pCompl = reinterpret_cast<FftwComplex*>(arraysC(iGridC2R).dataFirst());

    ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() << ATLOGVAR(iGridC2R) \
		   << ATLOGVAR(pReal) << ATLOGVAR(pCompl) \
		   << ATLOGVAR(*pReal) \
		   << ATLOGVAR((*pCompl)[0]) << ATLOGVAR((*pCompl)[1]) << "\n");

    fftwPlanC2R.dft_c2r(fftSize.length(),fftSize.dataFirst(),
			pCompl,
			pReal,
			FFTW_PATIENT);			    

    ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() << ATLOGVAR(sizeof(FftwComplex)) 
		   << ATLOGVAR(sizeof(arraysC(iGridC2R).data()[0])) << "\n");

  }

  Grid& getGrid(int ind) {

    return grids(ind);

  }

  void preprocessReceptor() {

    fftwPlansR2C(iGridRec).execute();

  }

  // The programm will spend most of its time here

  void correlate() {

    fftwPlansR2C(iGridLig).execute();
    arraysC(iGridLig) *= blitz::conj(arraysC(iGridRec));
    fftwPlanC2R.execute();
    //TODO:
    // Maybe optimize for speed by moving normalization to after the selection stage,
    // in other words, divide only the selected values
    int N = blitz::product(fftSize);
    grids(iGridC2R).getGridArray() /= N;

  }


  Grid
  unwrap(int indGrid) {

    typedef Math::WrappedIndex<T_num,N_dim> WrappedIndexType;

    dbg::trace t1(DBG_HERE);

    Grid& gridIn = getGrid(indGrid);

    WrappedIndexType wrappedIndex(fftSize);

    Grid gridUnwrapped(gridIn);

    gridUnwrapped.makeUnique();

    // wrappedIndex.unwrap() will center output at fftSize/2,
    // we make sure that gridUnwrapped will have spatial coord 0
    // at this index

    ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() \
		   << ATLOGVAR(gridUnwrapped.getGeometry().toSpatial(fftSize/2)) \
		   << ATLOGVAR(gridUnwrapped.getGeometry().toSpatial(fftSize/2 - 1)) \
		   << ATLOGVAR(gridUnwrapped.getGeometry().toSpatial(fftSize/2 + 1)) \
		   << ATLOGVAR((fftSize/2)) \
		   << ATLOGVAR(fftSize) \
		   << ATLOGVAR(gridUnwrapped.getGeometry().toLogical(Point(0))) \
		   << "\n");

    ATLOG_SWITCH_1(dbg::assertion(dbg::error, \
				  DBG_ASSERTION( \
						blitz::all(\
							   blitz::abs( \
								      gridUnwrapped.getGeometry().toLogical(Point(0)) - \
								      (fftSize/2) \
								      ) \
							   <= IntPoint(1))\
						)));
	
    gridUnwrapped = 0.;

    // ... Loop over all indices ...

    wrappedIndex.unwrap(gridIn.getGridArray(),gridUnwrapped.getGridArray());
      
    return gridUnwrapped;

  }


  // Test that FFT*FFT^-1 is identity operation:
  // does forward followed by reverse FFT of a ligand,
  // followed by normalization.
  // The result should be the original ligand grid.

  void testIdentity() {

    fftwPlansR2C(iGridLig).execute();	
    fftwPlanC2R.execute();
    int N = blitz::product(fftSize);
    grids(iGridC2R).getGridArray() /= N;	

  }

  FFTCorrelator() {}

  FFTCorrelator(const PointPair& boxDiag, T_num gridStep) {

    this->init(boxDiag,gridStep);

  }

  FFTCorrelator(const FFTCorrelator& x):

    grids(x.grids),
    fftwPlansR2C(x.fftwPlansR2C),
    fftwPlanC2R(x.fftwPlanC2R),
    fftSize(x.fftSize) {

    ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

    for(int i = 0; i < arraysC.length(); i++)
      arraysC(i).reference(x.arraysC(i));

  }

  FFTCorrelator& operator= (const FFTCorrelator& x) {

    ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

    if(this != &x) {

      //BUG:? copying grids like below is probably wrong
      //(elements should call reference())

      grids = x.grids;
      fftwPlansR2C = x.fftwPlansR2C;
      fftwPlanC2R = x.fftwPlanC2R;
      fftSize = x.fftSize;

      for(int i = 0; i < arraysC.length(); i++)
	arraysC(i).reference(x.arraysC(i));	  

    }
    return *this;
  }

  IntPoint sizeFft() const {

    return fftSize;

  }

protected:

  GridsSet grids;

  GridArraysCSet arraysC;

  FftwPlansSet fftwPlansR2C;
  FftwPlanType fftwPlanC2R;

  IntPoint fftSize;

};


typedef boost::shared_ptr<FFTCorrelator> PFFTCorrelator;

typedef typename common_types::num_vector_type<PFFTCorrelator>::Type VPFFTCorrelator;


class FFTCorrelators : public boost::noncopyable {

public:
      
  FFTCorrelators(const PointPair& boxDiag, T_num gridStep, int nFfts) {

    ATLOG_ASSERT_1(nFfts >= 1);

    ffts.resize(nFfts);

    gridsR.resize(nFfts);

    gridsL.resize(nFfts);

    gridsO.resize(nFfts);

    for(int i = 0; i < nFfts; i++) {
	  
      ffts(i).reset(new FFTCorrelator(boxDiag,gridStep));

      gridsR(i) = &(ffts(i)->getGrid(FFTCorrelator::iGridRec));
      gridsL(i) = &(ffts(i)->getGrid(FFTCorrelator::iGridLig));
      gridsO(i) = &(ffts(i)->getGrid(FFTCorrelator::iGridOut));

    }

    gridT = gridsO(0);

  }

  int size() const {

    return ffts.size();

  }

  void preprocessReceptor() {

    for( int i = 0; i < this->size(); i++ ) {
	  
      ffts(i)->preprocessReceptor();

    }
  }

  void correlate() {

    for( int i = 0; i < this->size(); i++ ) {
	  
      ffts(i)->correlate();

    }	
  }


  IntPoint sizeFft() const {

    return ffts(0)->sizeFft();

  }

  VPFFTCorrelator& getFfts() {

    return ffts;

  }


  VPRawGrid& getGridsRec() {

    return gridsR;

  }

  VPRawGrid& getGridsLig() {

    return gridsL;

  }

  VPRawGrid& getGridsOut() {

    return gridsO;

  }
      

  PRawGrid getGridTot() {

    return gridT;

  }


protected:

  // Shared ptrs to FFTs

  VPFFTCorrelator ffts;

  // Raw pointers to Receptor grids

  VPRawGrid gridsR;

  // Raw pointers to Ligand grids

  VPRawGrid gridsL;

  // Raw pointers to output grids (might be gridsL)

  VPRawGrid gridsO;

  // Raw pointer to total potential grid (might be one of gridsO)

  PRawGrid   gridT;

};


typedef boost::shared_ptr<FFTCorrelators> PFFTCorrelators;



#endif // PRODDL_DOCKING_FFT_H__

