//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_ROTATIONAL_GRID_H__
#define PRODDL_ROTATIONAL_GRID_H__

#include <fstream>
#include <string>

#include "PRODDL/Geom/transformation.hpp"
#include "PRODDL/Common/common_types.hpp"
#include "PRODDL/Common/debug.hpp"

#include "PRODDL/Common/logger.hpp"

#include "boost/filesystem/path.hpp"

namespace PRODDL { namespace Geom {

	template<typename T_num>
	typename TransformationTraits<T_num>::VRotation
		loadRotationalGrid(const std::string& filename) {

			ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

			typedef typename TransformationTraits<T_num>::VRotation VRot;
			typedef typename SpaceTraits<T_num>::Point3 Point;
			std::ifstream in(filename.c_str());
			ATALWAYS(in.good(),"loadRotationalGrid(): Unable to open input file: "+filename);
			int n_rotations = -1;
			in >> n_rotations;
			ATALWAYS(n_rotations > 0,"loadRotationalGrid(): Unable to read the number of rotations from file: " + filename);
			ATLOG_SWITCH_3(dbg::out(dbg::info) << dbg::indent() << "Creating " << n_rotations << " rotations...\n");
			VRot rotations(n_rotations);
			ATLOG_SWITCH_3(dbg::out(dbg::info) << dbg::indent() << ATLOGVAR(rotations.size()) << "\n");
			in.ignore(1000,'\n');
			T_num Pi = T_num(4.0)*std::atan(T_num(1.0));
			T_num degToRad = Pi/180.;

			ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() << ATLOGVAR(degToRad) << "\n");

			for(int i=0; i < n_rotations; i++) {

				Point angles;
				in >> angles(0) >> angles(1) >> angles(2);

				dbg::assertion(dbg::error, DBG_ASSERTION(in.good() || (i == (n_rotations - 1))));

				angles *= degToRad;

				if(i >= (n_rotations - 5)) {
					ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() << ATLOGVAR(i) << ATLOGVAR(angles) << "\n"); 
				}

				rotations(i) = Rotation<T_num>(angles);

			}
			return rotations;
	}

	template<typename T_num>
	typename TransformationTraits<T_num>::VRotation
		loadRotationalGrid(const std::string& dirname, int angle_step) {
			ATALWAYS(angle_step < 100 && angle_step > 0,"Angle step must be in (0,100) interval.");
			std::ostringstream out;
			out << "ang"
				<< std::setw(2) << std::setprecision(2) 
				<< std::setfill('0') << angle_step
				<< ".dat";
			std::string s_file = (boost::filesystem::path(dirname,boost::filesystem::native)/
				out.str()).PRODDL_BOOST_FILE_STRING();
			return loadRotationalGrid<T_num>(s_file);
	}

	//This version will look for the nearest existing value of rotational grid step

	template<typename T_num>
	typename TransformationTraits<T_num>::VRotation
		loadRotationalGrid(const std::string& dirname, int angle_step, int& angle_step_found) {

			ATTRACE_SWITCH_3(dbg::trace t1(DBG_HERE));

			//TODO: implement it in a sane manner - using boost:dir for instance.
			//Currently it tries to open files in a loop to find the nearest existing 
			//angle step file.
			if( ! (angle_step < 100 && angle_step > 0) ) {
				int req_step = angle_step;
				if(angle_step <=0)
					angle_step = 5;
				if(angle_step >= 100)
					angle_step = 90;
				ATLOG_SWITCH_1(dbg::out(dbg::warning) << dbg::indent() \
					<< "Requested angle step is " << req_step << ". " \
					<< "It must be in (0,100) degrees interval. Value changed to " \
					<< angle_step << "\n");
			}
			std::string rotGridFile;
			int angle_step_try = 0;
			bool next = true;
			for(int i_try = 0; i_try < 50 && next; i_try++) {
				for(int delta_sign = -1; delta_sign <= 1 && next; delta_sign += 2) {
					angle_step_try = angle_step + i_try*delta_sign;
					if(angle_step_try < 100 && angle_step_try > 0) {
						std::ostringstream out;
						out << "ang"
							<< std::setw(2) << std::setprecision(2) 
							<< std::setfill('0') << angle_step_try
							<< ".dat" << std::ends;
						std::string s_file = (boost::filesystem::path(dirname)/
							out.str()).string();

						ATLOG_SWITCH_3(dbg::out(dbg::info) << dbg::indent() \
							<< ATLOGVAR(angle_step_try) << ATLOGVAR(s_file) << "\n");
						std::ifstream in(s_file.c_str());
						if( in.good() ) {
							rotGridFile = s_file;
							next = false;
						}
					}
				}
			}
			ATLOG_SWITCH_1(dbg::out(dbg::info) << dbg::indent() << ATLOGVAR(rotGridFile) << "\n");
			ATALWAYS( ! rotGridFile.empty(), "loadRotationalGrid: No suitable rotational grid file located.");
			angle_step_found = angle_step_try;
			return loadRotationalGrid<T_num>(rotGridFile);
	}

}} // namespace PRODDL::Geom

#endif // PRODDL_ROTATIONAL_GRID_H__
