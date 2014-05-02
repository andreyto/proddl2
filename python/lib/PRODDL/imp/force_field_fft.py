"""Methods to load force field components needed for docking FFT scan"""
import IMP.atom
import IMP.container
import IMP.core

import logging, sys, datetime
log = logging.getLogger(__file__)


def extract_molforce(pdb_files):
    # Create an IMP model and add a heavy atom-only protein from a PDB file
    m = IMP.Model()

    # Read in the CHARMM heavy atom topology and parameter files
    ff = IMP.atom.get_heavy_atom_CHARMM_parameters()

    for i_prot,pdb_file in enumerate(pdb_files):
        prot = IMP.atom.read_pdb(pdb_file, m,
                                 IMP.atom.NonWaterNonHydrogenPDBSelector())


        # Using the CHARMM libraries, determine the ideal topology (atoms and their
        # connectivity) for the PDB file's primary sequence
        topology = ff.create_topology(prot)

        # Typically this modifies the C and N termini of each chain in the protein by
        # applying the CHARMM CTER and NTER patches. Patches can also be manually
        # applied at this point, e.g. to add disulfide bridges.
        topology.apply_default_patches()

        # Make the PDB file conform with the topology; i.e. if it contains extra
        # atoms that are not in the CHARMM topology file, remove them; if it is
        # missing atoms (e.g. sidechains, hydrogens) that are in the CHARMM topology,
        # add them and construct their Cartesian coordinates from internal coordinate
        # information.
        topology.setup_hierarchy(prot)


        # Add non-bonded interaction (in this case, Lennard-Jones). This needs to
        # know the radii and well depths for each atom, so add them from the forcefield
        # (they can also be assigned manually using the XYZR or LennardJones
        # decorators):
        ff.add_radii(prot)
        ff.add_well_depths(prot)

        # Get a list of all atoms in the protein, and put it in a container
        atoms = IMP.atom.get_by_type(prot, IMP.atom.ATOM_TYPE)

        for atom_p in atoms:
            atom = atom_p.get_as_atom()
            yield dict(
                 file_name=pdb_file,
                 mol_name=pdb_file,
                 eps=ff.get_epsilon(atom),
                 rad=ff.get_radius(atom),
                 xyz=atom.get_as_xyz().get_coordinates(),
                 mass=atom.get_as_mass().get_mass()
                 )

