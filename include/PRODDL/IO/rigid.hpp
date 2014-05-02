//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_IO_RIGID_H__
#define PRODDL_IO_RIGID_H__

// Top level classes for docking protocol

#include <cstdio>
#include <fstream>

#include "PRODDL/types.hpp"

#include "PRODDL/exceptions.hpp"

namespace PRODDL {


	template<typename T_num>
	class IORigid {

	public:

		enum { SIGNATURE = 2005052501 };

		typedef Types<T_num> TypesT;

		typedef typename TypesT::Point Point;

		typedef typename TypesT::Rotation Rotation;

		typedef typename TypesT::Translation Translation;

		typedef typename TypesT::RotationTranslation RotationTranslation;

		typedef typename TypesT::RotTranValue RotTranValue;


		// binary output will be done in this precision
		typedef float T_IO_num;

		//TODO: make this type platform dependent so that it is always 32 bit
		//(int on IA64 is 64 bit, for example)

		typedef int Int_IO;

		typedef Types<T_IO_num> TypesT_IO;

		typedef typename TypesT_IO::Point PointIO;

		struct IORecord {
			T_IO_num value;
			PointIO ang;
			PointIO xyz;
		};

		typedef typename common_types::num_vector_type<IORecord>::Type IORecords;

	protected:

		IORecords ioBuffer;
		Int_IO n;

	public:

		IORigid() : n(0) {
		}

		template<class InpIter>
		void
			putCoords(InpIter start, int _n) {

				n = Int_IO(_n);

				if (ioBuffer.rows() < n) {
					ioBuffer.resize(n);
				}

				for (Int_IO i = 0; i < n; i++, ++start) {

					const RotTranValue& x = *start;

					IORecord& r = ioBuffer(i);

					r.value = x.value;

					Point ang, xyz;

					x.tran.anglesAndDisplacement(ang,xyz);

					for(int i = 0; i < ang.length(); i++ ) {
						r.ang(i) = T_IO_num(ang(i));
						r.xyz(i) = T_IO_num(xyz(i));
					}

				}

		}


		template<class InpIter>
		void
			writeCoords(const std::string& fileName, InpIter start, int _n, char format = 'b') {
				putCoords(start,_n);
				writeCoords(fileName,format);
		}


		void
			writeCoords(const std::string& fileName, char format = 'b') {

				if( format == 'b' ) {

					std::FILE *f = std::fopen(fileName.c_str(), "wb");

					if( ! f ) {

						throw io_error("Cannot open file for writing: " + fileName);

					}

					Int_IO signature = SIGNATURE;

					std::fwrite((void*) &signature, sizeof(signature), 1, f);

					std::fwrite((void*) &n, sizeof(n), 1, f);

					int n_written = std::fwrite((void*)ioBuffer.dataFirst(),  sizeof(IORecord),  n,  f);

					if( n_written != n ) {

						std::fclose(f);

						throw io_error("Error when writing to file: " + fileName);

					}

					std::fclose(f);

				}

				else if( format == 't' ) {

					std::ofstream f(fileName.c_str());

					f << n << '\n';

					f << std::setprecision(10);
					f.setf(std::ios::showpoint);

					for( int i = 0; i < n; i++ ) {

						const IORecord& r = ioBuffer(i);

						f << r.value  << '\t'
							<< r.ang(0) << '\t' 
							<< r.ang(1) << '\t' 
							<< r.ang(2) << '\t'
							<< r.xyz(0) << '\t' 
							<< r.xyz(1) << '\t' 
							<< r.xyz(2) << '\n';

					}

				}

				else {

					throw argument_error("Unknown argument value for 'format' in 'IORigid::writeCoords': "+format);

				}

		}


		int
			readCoords(const std::string& fileName, char format = 'b') {

				if( format == 'b' ) {

					std::FILE *f = std::fopen(fileName.c_str(), "rb");

					if( ! f ) {

						throw io_error("Cannot open file for reading: " + fileName);

					}

					Int_IO signature = 0;

					std::fread((void*) &signature, sizeof(signature), 1, f);

					if ( signature != SIGNATURE ) {

						throw io_error("Format signature mismatch for binary input file in 'IORigid::readCoords()': "
							+ fileName);

					}

					std::fread((void*) &n, sizeof(n), 1, f);

					if( n < 0 ) {

						throw io_error("Loaded value for a number of records is less than zero: " + fileName);

					}

					if (ioBuffer.rows() < n) {
						ioBuffer.resize(n);
					}

					int n_read = std::fread((void*)ioBuffer.dataFirst(),  sizeof(IORecord),  n,  f);

					if( n_read != n ) {

						std::fclose(f);

						throw io_error("Error when reading in 'IORigid::readCoords' from file: " + fileName);

					}

					std::fclose(f);

				}

				else if( format == 't' ) {

					std::ifstream f(fileName.c_str());

					f >> n;

					if( n < 0 ) {

						throw io_error("Loaded value for a number of records is less than zero: " + fileName);

					}

					if (ioBuffer.rows() < n) {
						ioBuffer.resize(n);
					}

					for( int i = 0; i < n; i++ ) {

						IORecord& r = ioBuffer(i);

						f >> r.value
							>> r.ang(0) >> r.ang(1) >> r.ang(2)
							>> r.xyz(0) >> r.xyz(1) >> r.xyz(2);

						if( ! f ) {

							throw io_error("Error while reading in 'IORigid::readCoords' from file" + fileName);

						}

					}

				}

				else {

					throw argument_error("Unknown argument value for 'format' in 'IORigid::readCoords': "+format);

				}

				return int(n);

		}


		// This must be called after calling readCoords(). The code calling readCoords() can use
		// the returned number of records to prepare output range for getCoords()

		template<class OutIter>
		OutIter
			getCoords(int inp_start, int inp_n, OutIter out) {

				Int_IO inp_start_io = Int_IO(inp_start);
				Int_IO inp_end_io = Int_IO(inp_start+inp_n);

				if ( inp_end_io > ioBuffer.rows() ) {

					throw size_error("'ioBuffer' size is less than requested number of records in 'IORigid::getCoords()'");

				}

				for (Int_IO i = inp_start_io; i < inp_end_io; i++, ++out) {

					RotTranValue x;

					const IORecord& r = ioBuffer(i);

					x.value = r.value;

					x.tran = RotationTranslation(r.ang,r.xyz);

                    *out = x;

				}
                
                return out;

		}

		void convertCoords(const std::string& inpFile, const std::string outFile, char inpFormat = 'b') {

			readCoords(inpFile, inpFormat);

			char outFormat = ( inpFormat == 'b' ? 't' : 'b' );

			writeCoords(outFile, outFormat);

		}


	}; // class IORigid


} // namespace PRODDL

#endif // PRODDL_IO_RIGID_H__
