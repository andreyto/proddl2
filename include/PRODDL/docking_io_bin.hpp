//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_DOCKING_H__
#  error DOCKING_IO_BIN.HPP MUST BE INCLUDED FROM WITHIN DOCKING.HPP
#endif

#ifndef PRODDL_DOCKING_IO_BIN_H__
#define PRODDL_DOCKING_IO_BIN_H__


/// Binary dump IO module for docking data structures
/// This is to be used for data exchange withing a single parallel run
/// compiled for the same architecture

class RotFftScanIO_Bin {
public:
	enum { SIGNATURE = 2014011401 };
	RotFftScanIO_Bin(const std::string& file_name, std::fstream::openmode mode) {
		ATLOG_TRACE_3;
		m_io = new std::fstream(file_name.c_str(),mode | std::fstream::binary);
        ATALWAYS(m_io->good(),"Could not open input file: "+file_name);
	}

	~RotFftScanIO_Bin() {
		delete m_io;
	}

	bool read_record(Rotation& rot, TranValues& tran_vals) {
		ATLOG_TRACE_3;
		if (m_io->good() ) {
			int signature=0;
			m_io->read(reinterpret_cast<char*>(&signature), sizeof(signature));
			if( m_io->good() ) {
				ATALWAYS(signature==SIGNATURE,"Signature mismatch for binary format");
				//m_io->exceptions ( std::iostream::failbit | std::iostream::badbit );
				m_io->read(reinterpret_cast<char*>(&rot), sizeof(rot));
				std::size_t n=0;
				m_io->read(reinterpret_cast<char*>(&n), sizeof(n));
				tran_vals.resize(n);
				ATALWAYS(tran_vals.isStorageContiguous(),"Need output arrays with contiguous storage");
				m_io->read(reinterpret_cast<char*>(tran_vals.dataFirst()), n*sizeof(typename TranValues::T_numtype));
				if( m_io->good() ) 
					return true;
			}
		}
		return false;
	}

	bool write_record(const Rotation& rot, const TranValues& tran_vals) {
		ATLOG_TRACE_3;
		if (m_io->good()) {
			m_io->exceptions ( std::iostream::failbit | std::iostream::badbit );
			int signature = SIGNATURE;
			m_io->write(reinterpret_cast<const char*>(&signature), sizeof(signature));
			m_io->write(reinterpret_cast<const char*>(&rot), sizeof(rot));
			std::size_t n=tran_vals.size();
			m_io->write(reinterpret_cast<const char*>(&n), sizeof(n));
			ATALWAYS(tran_vals.isStorageContiguous(),"Need output arrays with contiguous storage");
			m_io->write(reinterpret_cast<const char*>(tran_vals.dataFirst()), n*sizeof(typename TranValues::T_numtype));
			return true;
		}
		else {
			return false;
		}
	}


protected:
	std::iostream *m_io;
};

// Reader to iterate over multiple instances of RotFftScanIO

template<class IO>
class RotFftScanIO_Collector {
public:

	RotFftScanIO_Collector(const std::string& file_name) {
		ATLOG_TRACE_3;
		m_io_list.reset(new std::fstream(file_name.c_str(),std::fstream::in));
		ATLOG_OUT_5("Opened list file: " << file_name);
		//m_io_list->exceptions ( std::iostream::failbit | std::iostream::badbit );
	}

	bool read_record(Rotation& rot, TranValues& tran_vals) {
		ATLOG_TRACE_3;
		bool status = false;
		if( m_io.get() ) {
		    ATLOG_OUT_5("Trying to read next fft scan record");
			status = m_io->read_record(rot,tran_vals);
		}
		if(! status) {
		    ATLOG_OUT_5("Could not read next fft scan record, checking next input file");
			if( m_io_list->good() ) {
                std::string file_name;
				if(std::getline(*m_io_list,file_name)) {
					m_io.reset(new IO(file_name,std::ios::in));
		            ATLOG_OUT_5("Opened new input file: " << file_name);
					status = m_io->read_record(rot,tran_vals);
				}
			}
		}
		ATLOG_OUT_5("Record status: " << status);
		return status;
	}

protected:
	boost::scoped_ptr<std::istream> m_io_list;
	boost::scoped_ptr<IO> m_io;
};



#endif // PRODDL_DOCKING_IO_BIN_H__

