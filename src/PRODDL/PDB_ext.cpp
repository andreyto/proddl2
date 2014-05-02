//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "PRODDL/Os/compressed_io.hpp"

#include <fstream>

#include "PRODDL/External/Pdb++/pdb++.hpp"

#include "PRODDL/Common/logger.hpp"

#include "PRODDL/Common/to_string.hpp"

#include "PRODDL/Geom/traits.hpp"

#include <map>

#include <vector>

namespace PRODDL {

  // Save one chain from PDB Biounit file as a separate file.
  // (PDB Biounit files use the same schema as PDB NMR files -
  // multiple chains in multiple models).
  // Return 1 if any ATOM records were written, 0 - otherwise.

  int extractBiounitChain(const std::string& inpFile,
			  const char chainId,
			  int modelId,
			  const std::string& outFile) {

    typedef PDBPP::PDB Record;
    std::auto_ptr<std::istream> p_in = 
      openUncompressedInput(inpFile.c_str());

    ATLOG_SWITCH_1(dbg::assertion(dbg::error, \
				  DBG_ASSERTION( p_in->good() )));

    std::ofstream out(outFile.c_str());

    ATLOG_SWITCH_1(dbg::assertion(dbg::error, \
				  DBG_ASSERTION( out.good() )));

    Record record;

    int curModelId = 0;

    bool copy_record = false;

    int n_found = 0;

    while ( (*p_in) >> record ) {

      if ( record.type() == PDBPP::PDB::MODEL ) {
	curModelId = record.model.num;
	copy_record = true;
      }
      else if ( record.type() == PDBPP::PDB::ENDMDL ) {
	curModelId = -1;
	copy_record = true;
      }
      else if ( record.type() == PDBPP::PDB::END ) {
	curModelId = -1;
	copy_record = true;
      }
      else if ( record.type() == PDBPP::PDB::ATOM &&
		record.atom.residue.chainId == chainId && 
		curModelId == modelId ) {

	copy_record = true;
	n_found = 1;
      }
      else {
	copy_record = false;
      }
      if ( copy_record ) {
	out << record << '\n';
      }
    }
    return n_found;
  }


  // Save certain chains from certain models from PDB Biounit file as a separate file.
  // Will select chains identified by (chainIds[i],modelIds[i]).
  // (PDB Biounit files use the same schema as PDB NMR files -
  // multiple chains in multiple models).
  // Return 1 if any ATOM records were written, 0 - otherwise.

  int extractBiounitChains(const std::string& inpFile,
			   boost::python::dictionary selection,
			   const std::string& outFile) {

    namespace python = boost::python;

    typedef std::map<std::string,std::string> Smap;

    Smap sel;

    python::list items = selection.items();

    for( int i_item = 0; i_item < items.size(); i_item++) {

      python::tuple item(items.get_item(i_item));

      std::string key = python::string(item[0]).c_str();
      std::string val = python::string(item[0]).c_str();

      sel[key] = val;
    }

    typedef PDBPP::PDB Record;
    std::auto_ptr<std::istream> p_in = 
      openUncompressedInput(inpFile.c_str());

    ATLOG_SWITCH_1(dbg::assertion(dbg::error, \
				  DBG_ASSERTION( p_in->good() )));

    std::ofstream out(outFile.c_str());

    ATLOG_SWITCH_1(dbg::assertion(dbg::error, \
				  DBG_ASSERTION( out.good() )));

    Record record;

    int curModelId = 0;

    bool copy_record = false;

    int n_found = 0;

    while ( (*p_in) >> record ) {

      copy_record = false;

      if ( record.type() == PDBPP::PDB::MODEL ) {
	curModelId = record.model.num;
      }
      else if ( record.type() == PDBPP::PDB::ENDMDL ) {
	curModelId = -1;
      }
      else if ( record.type() == PDBPP::PDB::END ) {
	curModelId = -1;
	copy_record = true;
      }
      else if ( record.type() == PDBPP::PDB::ATOM || 
		record.type() == PDBPP::PDB::TER ) {
	Record::Residue *p_res = 0;
	if ( record.type() == PDBPP::PDB::ATOM ) {
	  p_res = &record.atom.residue;
	}
	else {
	  p_res = &record.ter.residue;
	}
	std::string key((p_res->chainId + to_string(curModelId)).c_str());
	if( sel.find(key) != sel.end() ) {
	  copy_record = true;
	  n_found = 1;
	  p_res->chainId = sel[key][0];
	}
      }

      if ( copy_record ) {
	out << record << '\n';
      }
    }
    return n_found;
  }


  template<typename T_num>
  class PDBConfiguration {

  public:


    typedef PDBPP::PDB Record;

    typedef std::vector<Record> Records;

    typedef typename Geom::SpaceTraits<T_num>::Point3 Point;

    typedef typename Geom::SpaceTraits<T_num>::VPoint3 Points;

    typedef common_types::num_index_type::Type Index;


    PDBConfiguration() {}

    PDBConfiguration(const std::string& fileName) {

      //ATTRACE_SWITCH_4(dbg::trace t1(DBG_HERE));

      load(fileName);

    }

    void load(const std::string& fileName) {
      //TODO: add flag that would cause the validation of structure
      //identity between the old and newly loaded files:
      //by ATOM record comparison - name and relative index must be the same

      int sizeReserve = records.size();

      if ( sizeReserve < 1000 )
	sizeReserve = 1000;

      records.clear();

      ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() \
		     << ATLOGVAR(sizeReserve) << '\n');

      records.reserve(sizeReserve);

      std::auto_ptr<std::istream> p_in = openUncompressedInput(fileName.c_str());

      //std::ifstream *p_in = new std::ifstream(fileName.c_str());

      ATLOG_SWITCH_1(dbg::assertion(dbg::error, DBG_ASSERTION( p_in->good() )));

      Record record;

      // We drop all records except ATOM records. (SCWRL, for instance, can add SSBOND records
      // for different conformations).

      while ( (*p_in) >> record ) {

	bool skip = true;

	if ( record.type() == PDBPP::PDB::ATOM ) {

	  // for now, we just skip the alternative locations

	  if ( record.atom.altLoc == ' ' )
	    skip = false;

	  ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() << ATLOGVAR(record.atom.altLoc) \
			 << ATLOGVAR(record.atom.xyz[2]) << '\n');

	}
      
	ATLOG_SWITCH_4(dbg::out(dbg::info) << dbg::indent() \
		       << ATLOGVAR(skip) << ATLOGVAR(record) << '\n');

	if ( ! skip )
	  records.push_back(record);

      }

    }


    void save(const std::string& fileName) {

      std::auto_ptr<std::ostream> p_out = openCompressedOutput(fileName.c_str());

      ATLOG_SWITCH_1(dbg::assertion(dbg::error, DBG_ASSERTION( p_out->good() )));

      for(int i = 0; i < records.size(); i++ ) {

	(*p_out) << records[i] << '\n';

      }

    }


    int countAtoms() const {

      int n_atoms = 0;

      for(int i = 0; i < records.size(); i++ ) {

	if( records[i].type() == PDBPP::PDB::ATOM )
	  n_atoms++;

      }

      return n_atoms;

    }


  protected:

    Records records;

  };


} // namespace PRODDL


namespace PRODDL {

  Logger gLogger;

} // namespace PRODDL


