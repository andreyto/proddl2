#include "PRODDL/External/Pdb++/pdb++.hpp"

#include <cstdlib>

#include <iostream>

//#include <gzstream.h>

#include <fstream>

#include <string>

int main() {

  std::cout << "Reading pdb file names from standard input." << std::endl;

  int n_files = 0;

  std::string file_name;

  while(std::cin >> file_name) {

    std::ifstream file_in(file_name.c_str());

    if( ! file_in ) {
      std::cerr << "Unable to open input pdb file: " << file_name << std::endl;
      continue;
    }

    PDBPP::PDB	record;

    int n_atoms = 0;

    while (file_in >> record) {
	std::cout << record << '\t' << record.chars() << "\n";
      switch (record.type()) {
      case PDBPP::PDB::ATOM:
	      std::cout << record.atom.altLoc << ' '
	                << record.atom.xyz[0] << ' ' << record.atom.xyz[1]
		        << ' ' << record.atom.xyz[2] << '\n';
	n_atoms++;
	break;
      case PDBPP::PDB::END:
	n_files++;
	break;
      }
    }

    std::cout << "Number of ATOM records: " << n_atoms << std::endl;

  }
 
  std::cout << "Number of END  records: " << n_files << std::endl;

  return 0;
}
