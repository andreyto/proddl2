//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
// Handling of files in Protein Data Bank format

#ifndef PRODDL_PDB_H__
#define PRODDL_PDB_H__

#include <pdb++.hpp>

namespace PRODDL {

  class Pdb
  {
  public:

    // Load PDB file in memory, either fully or everything except the ATOM and HETATM records
    
    void loadFile(const std::string& fileName, bool load_coords);
    
  }; // class Pdb 

} // namespace PRODDL

#endif // PRODDL_PDB_H__
