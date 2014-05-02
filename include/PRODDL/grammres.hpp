//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef AT__GRAMMRES_H__
#define AT__GRAMMRES_H__

#include <vector>
#include <string>
#include <iostream>


#include "PRODDL/Geom/traits.hpp"

#include "PRODDL/Common/logger.hpp"

namespace PRODDL {

  // Read and keep gramm *.res file

  class gr_rmol_chain
  {
  public:
    std::string sid; //ID assigned to that subsystem
    std::string sfile; //pdb file
    std::string schain; //chain id
    int cat; //atom count
    gr_rmol_chain() : cat(0) {}
    void print(std::ostream& out) const
    {
      out << '\t' << sid << '\t' << sfile << '\t' 
	  << schain << '\t' << cat << std::endl;
    }
  };

  inline std::ostream& operator<< (std::ostream& out,const gr_rmol_chain& x)
  {
    x.print(out);
    return out;
  }

  class gr_match
  {
  public:

    // numeric type
    typedef double mfl;
    // 3D vector type
    typedef Geom::SpaceTraits<mfl>::Point3 vect;
    int ind; //index of the match (0 offset)
    mfl e; //"energy"
    vect ang; //euler angles
    vect xyz; //CM coords
    gr_match(): ind(0), e(0), ang(0,0,0), xyz(0,0,0) {}
    void print(std::ostream& out) const;
  };

  inline std::ostream& operator<< (std::ostream& out,const gr_match& x)
  {
    x.print(out);
    return out;
  }

  class gr_res
  {
  public:
    // type for array of matches
    typedef std::vector<gr_match> amatch;
    // type for pair of chains
    typedef std::pair<gr_rmol_chain,gr_rmol_chain> molpair;
    // numeric type
    typedef gr_match::mfl mfl;
    // 3D point type
    typedef gr_match::vect vect;

    amatch matches; //records of matches
    molpair chains; //what was docked
    mfl grdstep; //grid step
    int cmatch; //number of matches in output as in a field of res file
    std::string sfres; //name of res file
    vect initial_angles; // initial rotation angles (from [initial_angles] 
    //section, will be set to (0,0,0) if this section is 
    //not present in .res file. Upon reading matches from .res file,
    //they will be automatically converted to reflect the [initial_angles]
    //parameter, so no user processing of this parameter is normally
    //required, it is provided only for information
    gr_res() : grdstep(0), cmatch(0) {}
    explicit gr_res(const std::string& sfile) //construct from file
    {
      init(sfile);
    }
    explicit gr_res(const std::string& sfile, int max_matches) //construct from file, load no more than 'max_matches'
    {
      init(sfile, max_matches);
    }
    void init(const std::string& sfile, int max_matches = -1); //init from file
    void print(std::ostream& out) const ; //print out for debugging

    mfl get_grdstep() const { return grdstep; }
    const std::string& get_sfres() const { return sfres; }
    molpair& get_chains() { return chains; }
    const molpair& get_chains() const { return chains; }
    const amatch& get_matches() const { return matches; }
  };

  inline std::ostream& operator<< (std::ostream& out,const gr_res& x)
  {
    x.print(out);
    return out;
  }


} // namespace PRODDL

#endif // AT__GRAMMRES_H__
