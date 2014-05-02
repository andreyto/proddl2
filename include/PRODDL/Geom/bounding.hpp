//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_GEOM_BOUNDING_H__
#define PRODDL_GEOM_BOUNDING_H__

#include "PRODDL/Geom/transformation.hpp"
#include "PRODDL/Geom/traits.hpp"
#include "PRODDL/Geom/gdiam_simple.hpp"

#include "PRODDL/Common/common_types.hpp"
#include "PRODDL/Common/bz_cast.hpp"
#include "PRODDL/Blitz/bzutil.hpp"

#include "PRODDL/Blitz/bzalgor.hpp"

#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/logger.hpp"

#include "PRODDL/Common/common_algor.hpp"

namespace PRODDL { namespace Geom {

	namespace Bounding {

		template<typename T_num>
		class Diameter {

		public:

			typedef typename SpaceTraits<T_num>::Point3 Point;
			typedef typename SpaceTraits<T_num>::Point3Pair PointPair;
			typedef typename SpaceTraits<T_num>::VPoint3 VPoint;

		protected:

			typedef double T_num_internal;
			typedef typename SpaceTraits<T_num_internal>::Point3 PointInternal;
			typedef typename SpaceTraits<T_num_internal>::Point3Pair PointPairInternal;
			typedef typename SpaceTraits<T_num_internal>::VPoint3 VPointInternal;
			typedef typename SpaceTraits<T_num_internal>::LargeMatrix LargeMatrixInternal;

			PointPair diam;

		public:

			Diameter() {}

			Diameter(const VPoint& _points, T_num eps=0.0) {
				ATLOG_TRACE_3;
				VPointInternal points(blitz_ext::getContiguousRefOrCopy<PointInternal>(_points));
				LargeMatrixInternal pointsAsMatrix(blitz::viewWithFlattenedComponent(points));
				PointPairInternal diamIntern;
				PRODDL::Geom::Gdiam::simp_approx_diam_pair(pointsAsMatrix.dataFirst(),
					pointsAsMatrix.rows(),
					diamIntern(0).data(),
					diamIntern(1).data(),
					T_num_internal(eps));

				//copy the function's output to T_num array
				for(int i=0; i < diam.length(); i++) {
					for(int j=0; j < diam(i).length(); j++) {
						diam(i)(j) = T_num(diamIntern(i)(j));
					}
				}

				ATLOG_OUT_3(ATLOGVAR(diam) << ATLOGVAR(diamIntern));

				//diam = blitz::cast(diamIntern,T_num());
			}

			const PointPair& getDiameter() const {
				ATLOG_TRACE_3;
				return diam;
			}

			const T_num getSize() const {
				ATLOG_TRACE_3;
				return T_num(Math::normSelf(diam[1] - diam[0]));
			}

			// Return Translation transformation object to the coordinate
			// system with axes parallel to the current system, and
			// a center defined by the center of this diameter.

			Translation<T_num>
				getTransformationToBoundCoordinates() const {
					ATLOG_TRACE_3;
					return Translation<T_num>(Point( - (diam(1) + diam(0))/2) );
			}

		};


		template<typename T_num>
		class Box {

		public:

			typedef typename SpaceTraits<T_num>::Point3 Point;
			typedef typename SpaceTraits<T_num>::Point3Pair PointPair;
			typedef typename SpaceTraits<T_num>::Matrix3x3 Matrix3;
			typedef typename SpaceTraits<T_num>::VPoint3 VPoint;


		protected:

			typedef double T_num_internal;
			typedef typename SpaceTraits<T_num_internal>::Point3 PointInternal;
			typedef typename SpaceTraits<T_num_internal>::Point3Pair PointPairInternal;
			typedef typename SpaceTraits<T_num_internal>::Matrix3x3 Matrix3Internal;
			typedef typename SpaceTraits<T_num_internal>::VPoint3 VPointInternal;
			typedef typename SpaceTraits<T_num_internal>::LargeMatrix LargeMatrixInternal;

			// main diagonal of the box
			PointPair diag;

			// directions of box sides, in rows of the matrix
			Matrix3 directions;

		public:

			Box() {}

			// 'fixedDirections' - if true, the box will have orientation along the current
			// coordinate axes - the diagonal will be essentially the min and max
			// coordinates in the point set;
			// Otherwise, the non-trivial operation will be performed, which will construct
			// the smallest possible box in the best orientation.

			Box(const VPoint& _points, bool fixedDirections=false,int grid_n=5, int _sample_n=-1) {
				ATLOG_TRACE_3;
				if( ! fixedDirections ) {
					VPointInternal points(blitz_ext::getContiguousRefOrCopy<PointInternal>(_points));

					//ATLOG_OUT_4(ATLOGVAR(_points) << ATLOGVAR(points));

					LargeMatrixInternal pointsAsMatrixIntern(blitz::viewWithFlattenedComponent(points));
					PointPairInternal diagIntern;
					Matrix3Internal directionsIntern;
					// negative value of _sample_n means we want to sample all points
					int sample_n = _sample_n < 0 ? pointsAsMatrixIntern.rows() : _sample_n;
					PRODDL::Geom::Gdiam::simp_approx_mvbb_grid_sample(pointsAsMatrixIntern.dataFirst(),
						pointsAsMatrixIntern.rows(),
						diagIntern(0).data(),
						diagIntern(1).data(),
						directionsIntern.data(),
						grid_n, sample_n);
					//copy the function's output to T_num arrays
					for(int i=0; i < diag.length(); i++) {
						for(int j=0; j < diag(i).length(); j++) {
							diag(i)(j) = T_num(diagIntern(i)(j));
						}
					}

					ATLOG_OUT_3(ATLOGVAR(diag) << ATLOGVAR(diagIntern));

					//diag = blitz::cast(diagIntern,T_num());

					for(int i=0; i < 3; i++) {
						for(int j=0; j < 3; j++) {
							directions(i,j) = T_num(directionsIntern(i,j));
						}
					}

					ATLOG_OUT_3(ATLOGVAR(directions) << ATLOGVAR(directionsIntern));

					//directions = blitz::cast(directionsIntern,T_num());
				}
				else {
					diag(0) = 0; diag(1) = 0;
					blitz_ext::bracketer _bracketer;
					int n_points = _points.size();
					if(n_points > 0) {
						diag(0) = _points(0);
						diag(1) = diag(0);
						for(int i = 0; i < n_points; i++) {
							_bracketer(_points(i),diag(0),diag(1));
						}
					}

					directions = 
						1,0,0,
						0,1,0,
						0,0,1;
				}
			}


			// Return RotationTranslation transformation object to the coordinate
			// system defined by the directions and center of this bounding box.

			RotationTranslation<T_num>
				getTransformationToBoundCoordinates() const {
					ATLOG_TRACE_3;
					return Rotation<T_num>(directions) *
						Translation<T_num>(Point( - (diag(0) + (diag(1) - diag(0))/2) ));
			}

			// Move this object to the coordinate system defined by its directions and center,
			// by applying transformation returned by getTransformationToBoundCoordinates() method.
			// Return the transformation.
			// After this method is called, the getDiagonal() method will return the pair of lower and upper
			// coordinates of the box, and the size of the box will be diagonal[1] - diagonal[0].

			RotationTranslation<T_num> moveToBoundCoordinates() {
				ATLOG_TRACE_3;
				RotationTranslation<T_num> transform = getTransformationToBoundCoordinates();
				applyTransformation(transform);
				ATALWAYS(blitz::all((diag(1) - diag(0)) >= 0),\
					"tranformationToBoxCoordinates(): Main diagonal in box coordinates must point strictly from low to high.");
				return transform;
			}

			// Move this object by applying coordinate transformation

			void applyTransformation(const RotationTranslation<T_num>& transform) {
				ATLOG_TRACE_3;
				diag(0) = transform(diag(0));
				diag(1) = transform(diag(1));
				directions = blitz_ext::product(directions,PRODDL::Math::transpose(transform.getTensor()));
			}

			const Matrix3& getDirections() const {

				ATLOG_TRACE_3;

				ATLOG_OUT_3(ATLOGVAR(directions));

				return directions;
			}

			const PointPair& getDiagonal() const {

				ATLOG_TRACE_3;

				ATLOG_OUT_3(ATLOGVAR(diag));

				return diag;
			}

		};


		// bounding box in a fixed coordinate frame
		// orientation of the box is the same as orientation of the coordinate system axes


		template<typename T_num>
		class BoxFixed {

		public:

			typedef typename SpaceTraits<T_num>::Point3Pair PointPair;

		protected:

			PointPair diag;

		public:

			BoxFixed() {}

		};

	} //namespace Bounding

}} // namespace PRODDL::Geom

#endif // PRODDL_GEOM_BOUNDING_H__
