//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_DOCKING_H__
#define PRODDL_DOCKING_H__

// Top level classes for docking protocol

#include <cmath>

#include "PRODDL/types.hpp"

#include "PRODDL/Math/fft.hpp"

#include "PRODDL/Math/fftw.hpp"

#include "PRODDL/Common/queue.hpp"

#include "PRODDL/Common/math.hpp"

#include "PRODDL/Common/bz_cast.hpp"

#include "PRODDL/Geom/bounding.hpp"

#include "PRODDL/Geom/symmetry.hpp"

#include "PRODDL/Geom/rotational_grid.hpp"

#include "PRODDL/Geom/cluster.hpp"

#include "PRODDL/potentials.hpp"

#include "PRODDL/Common/g_options.hpp"

#include "PRODDL/Common/logger.hpp"

#include "PRODDL/Common/common_algor.hpp"

#include "PRODDL/IO/hdf5.hpp"

#include "PRODDL/IO/rigid.hpp"

#include <deque>
#include <iterator>
#include <algorithm>
#include <cstddef>
#include <boost/scoped_ptr.hpp>

namespace PRODDL {


	// Docking implements a concept of "parametrized namespace" - it
	// houses other classes which use its typedefs based on its
	// template arguments.
	// The included code is split into several '#include' files
	// for editing convenience. Thus, those '#include' files should
	// not themselves use '#include' constructs for standard headers - that must
	// be done at the top of this file. Also, those files must remember
	// that they already exist inside PRODDL::Docking<T_num> scope.


	template<typename T_num>
	class Docking {

	public:

		enum { N_dim = 3 };

		typedef Types<T_num> TypesT;

		typedef typename TypesT::Grid Grid;

		typedef typename TypesT::GridC GridC;

		typedef Grid* PRawGrid;

		typedef typename common_types::num_vector_type<PRawGrid>::Type VPRawGrid;

		typedef typename Grid::GridArray GridArray;

		typedef typename GridC::GridArray GridArrayC;

		typedef typename TypesT::Projector Projector;

		typedef typename TypesT::IntPoint IntPoint;

		typedef typename TypesT::Point Point;

		typedef typename TypesT::Points Points;

		typedef typename TypesT::Ints Ints;

		typedef typename TypesT::Floats Floats;

		typedef typename TypesT::Point3Pair PointPair;

		typedef typename TypesT::Rotation Rotation;

		typedef typename TypesT::Rotations Rotations;

		typedef std::deque<Rotation> DequeRotations;

		typedef typename TypesT::Translation Translation;

		typedef typename TypesT::Translations Translations;

		typedef typename TypesT::RotationTranslation RotationTranslation;

		typedef typename TypesT::RotationTranslations RotationTranslations;

		typedef typename TypesT::TranValue TranValue;

		typedef typename TypesT::TranValues TranValues;

		typedef typename TypesT::RotTranValue RotTranValue;

		typedef typename TypesT::RotTranValues RotTranValues;


		// coordinates of receptor and ligand

		// N_mol - number of molecules in docking (currently two body docking - receptor and ligand)
		// iRec and iLig - indices of receptor and ligand in various arrays

		enum { N_mol = 2, iRec = 0, iLig = 1 };

		// type to carry coordinates of molecules as one object

		typedef typename common_types::point_type<Points,N_mol>::Type PointsMolSet;

		// type to carry transformations of molecules as one object

		typedef typename common_types::point_type<RotationTranslation,N_mol>::Type TransformsMolSet;

		// type to carry floating point scalar arrays ( one for each molecule ) as one object

		typedef typename common_types::point_type<Floats,N_mol>::Type FloatsMolSet;

		typedef Potentials<T_num> PotentialsT;

		typedef typename PotentialsT::SoftCoreLJ LjPot;
		typedef typename PotentialsT::SoftCoreLJRep LjPotRep;
		typedef typename PotentialsT::SoftCoreLJAttr LjPotAttr;



#    include "PRODDL/docking_molforce.hpp"

#    include "PRODDL/docking_molstruct.hpp"

#    include "PRODDL/docking_fft.hpp"

#    include "PRODDL/docking_io_bin.hpp"

#    include "PRODDL/docking_app.hpp"

#    include "PRODDL/docking_io_hdf5.hpp"



	}; // class Docking


} // namespace PRODDL

#endif // PRODDL_DOCKING_H__
