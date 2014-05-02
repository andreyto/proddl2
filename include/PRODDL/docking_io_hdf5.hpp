//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_DOCKING_H__
#  error DOCKING_IO_HDF5.HPP MUST BE INCLUDED FROM WITHIN DOCKING.HPP
#endif

#ifndef PRODDL_DOCKING_IO_HDF5_H__
#define PRODDL_DOCKING_IO_HDF5_H__


/// HDF5 IO module for docking data structures

static
void load_from_hdf5(const std::string& file_name, MolForceParams& self) {
	ATLOG_TRACE_3;
	using namespace Hdf5;
	using namespace blitz;
	HDF5File::mode_t flag = HDF5File::in;
	HDF5File inp(file_name, flag);      
	inp.getAttribute("/","mix",self.mix);
	Array<int,2> mol_offsets(inp.readArray<int,2>("mol_offsets"));
	ATALWAYS(mol_offsets.rows() >= N_mol && mol_offsets.columns() == 2,"mol_offsets must be Nx2 matrix");
	Array<T_num,2> pos(inp.readArray<T_num,2>("pos"));
	Array<T_num,1> mass(inp.readArray<T_num,1>("mass"));
	Array<T_num,1> eps(inp.readArray<T_num,1>("eps"));
	Array<T_num,1> sigma(inp.readArray<T_num,1>("sigma"));
	Array<T_num,1> alpha(inp.readArray<T_num,1>("alpha"));

	for(int i_mol = 0; i_mol < N_mol; i_mol++) {

		int start = mol_offsets(i_mol,0);
		int end = mol_offsets(i_mol,1);
		Range range(start,end-1); //Blitz::Range is a closed range

		ATALWAYS(end <= alpha.ubound(0) + 1 && start >= alpha.lbound(0),\
			"Molecule offset is past array boundaries");

		self.pos(i_mol).reference( 
			viewWithFoldedComponent<TinyVector<T_num,N_dim> >
			(pos(range,Range::all()))
			.copy());
		self.mass(i_mol).reference(mass(range));
		self.eps(i_mol).reference(eps(range));
		self.sigma(i_mol).reference(sigma(range));
		self.alpha(i_mol).reference(alpha(range));

		int n = self.pos(i_mol).size();

		ATALWAYS(self.mass(i_mol).size() == n && 
			self.mass(i_mol).size() == n &&
			self.eps(i_mol).size()  == n &&
			self.sigma(i_mol).size()== n &&
			self.alpha(i_mol).size()== n,
			"Atom parameter arrays have different lengths");

		ATLOG_OUT_3(ATLOGVAR(self.pos(i_mol).size()) \
			<< ATLOGVAR(self.mass(i_mol).size()) \
			<< ATLOGVAR(self.eps(i_mol).size()) \
			<< ATLOGVAR(self.sigma(i_mol).size()) \
			<< ATLOGVAR(self.alpha(i_mol).size()) \
			<< ATLOGVAR(i_mol));

	}

}

#endif // PRODDL_DOCKING_IO_HDF5_H__

