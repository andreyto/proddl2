//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/grammres.hpp"

#include "PRODDL/Common/to_string.hpp"

#include <functional>
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>

#include <math.h>
#include <stdio.h>

#include "PRODDL/Common/debug.hpp"
#include "PRODDL/Common/string_util.hpp"
#include "PRODDL/Common/math.hpp"
#include "PRODDL/Geom/transformation.hpp"

namespace PRODDL {

  void gr_match::print(std::ostream& out) const
  {
    out << '\t' << ind << '\t' << e << '\t' 
	<< ang << '\t' << xyz << std::endl;
  }

  //Service routine for gr_res constructor.
  //Read one entry of [molecules] section fromstring s,
  //into gr_rmol_chain object ch.

  void parseChain(std::istream& in, gr_rmol_chain& ch)
  {
    //Sample input:
    //2ptcE  ( pdb2ptc.ent fragment E  1629 atoms)

    std::string ignore;

    in >> ch.sid >> ignore >> 
      ch.sfile >> ignore >> 
      ch.schain >> 
      ch.cat >> ignore;
  }

  void gr_res::init(const std::string& sfile,int max_matches)
  {
    ATLOG_TRACE_4;

    sfres = sfile;
    std::ifstream in(sfile.c_str());
    ATLOG_SWITCH_4(dbg::out(dbg::info) << \
		   dbg::indent() << ATLOGVAR(sfile) << "\n");
    ATLOG_SWITCH_1(dbg::assertion(dbg::error,\
				  DBG_ASSERTION( in.good() )));

    matches.clear();

    for(std::string word; in.good() && in >> word; ) {

      if( word == "[molecules]" ) {
	
	parseChain(in,chains.first);
	parseChain(in,chains.second);

	ATLOG_SWITCH_1(dbg::assertion(dbg::error,\
				      DBG_ASSERTION( in.good() )));

      }

      else if( word == "[initial_angles]" ) {

	in \
	  >> initial_angles[0]
	  >> initial_angles[1]
	  >> initial_angles[2];

	ATLOG_SWITCH_1(dbg::assertion(dbg::error,\
				      DBG_ASSERTION( in.good() )));

      }

      else if( word == "[match]" ) {

	for(std::string line; in.good() && std::getline(in,line); )
	  {
	    if( max_matches >= 0 && matches.size() >= max_matches )
	      break;
	    if(!strempty(line)) {
	      
	      gr_match m;
	      std::istringstream isl(line);
	      isl >> m.ind >> m.e 
		  >> m.ang[0] >> m.ang[1] >> m.ang[2] 
		  >> m.xyz[0] >> m.xyz[1] >> m.xyz[2];
	      if(m.ind > 0) { //did read something
		  
		//bring index to zero offset
		m.ind--; 
		//bring angles to radian units
		m.ang *= ( M_PI/180.0 ); 
		// distances to nanometers
		//m.xyz *= Math::Constants<mfl>::AngstromToNm; 
		//Here we transparently handle 
		//the [initial_angles] parameter:
		//convert m.ang to Rotation object, 
		//compose it with Rotation defined by 
		//[initial_angles], convert the
		//composition back to Euler angles and 
		//store them into m.ang
		typedef PRODDL::Geom::Rotation<mfl> Rotation;
		Rotation match_rot(m.ang);
		Rotation initial_rot(initial_angles); 
		Rotation final_rot = match_rot*initial_rot;
		//DEBUG:
		vect old_ang = m.ang;
		m.ang = final_rot.eulerAngles();
		//DEBUG:
// 		if( m.ind < 10) {
// 		  std::cout << "gramm.res match: " << m.ind << 
// 		    '\t' << old_ang << '\t' << m.ang 
// 			    << '\t' << m.xyz << std::endl;
// 		  //std::cout << Rotation(old_ang).getTensor() 
// 		  //<< std::endl << Rotation(m.ang).getTensor() 
// 		  //<< std::endl;
// 		} // if ( m.ind < 10
		matches.push_back(m);
	      } // if ( m.ind > 0
	    } // if ( ! strempty(line
	  }  // for ( line
      } // if ( word
    } // for ( word

    grdstep = 0;

    cmatch = matches.size();

  }


  void gr_res::print(std::ostream& out) const //print out for debugging
  {
    out << "Molecules:" << std::endl << chains.first 
	<< std::endl << chains.second << std::endl;
    out << "Grid step = " << grdstep << '\t'
	<< "Number of matches = " << cmatch << '\t'
	<< "initial_angles = " << initial_angles << '\t'
	<< std::endl;
    out << "Matches:" << std::endl; 
    std::copy(matches.begin(),matches.end(),
	      std::ostream_iterator<gr_match>(out,"\n"));
    out << std::endl;
  }

} // namespace PRODDL
