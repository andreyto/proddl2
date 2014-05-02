## \example of locally optimizing LJ potential between two rigid body proteins loaded from PDB
##
import IMP.atom
import IMP.container
import IMP.core
import IMP.gsl

import logging, sys, datetime
log = logging.getLogger(__file__)


class _log_formatter(logging.Formatter):
    converter=datetime.datetime.fromtimestamp
    def formatTime(self, record, datefmt=None):
        ct = self.converter(record.created)
        if datefmt:
            s = ct.strftime(datefmt)
        else:
            t = ct.strftime("%Y-%m-%d %H:%M:%S")
            s = "%s.%03d" % (t, record.msecs)
        return s

log_h_console = logging.StreamHandler()
log_h_console.setFormatter(_log_formatter(fmt='%(asctime)s %(message)s',datefmt='%Y-%m-%d,%H:%M:%S.%f'))
logging.root.addHandler(log_h_console)

#logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

# it gets awfully slow with internal checks
IMP.base.set_check_level(IMP.base.USAGE)
#IMP.base.set_log_level(IMP.base.SILENT)
IMP.base.set_log_timer(True)
lc = IMP.base.CreateLogContext(str(__file__+".log"))
lt = IMP.base.SetLogTarget(str(__file__+".log"))

#IMP.base.set_number_of_threads(4)

# Create an IMP model and add a heavy atom-only protein from a PDB file
m = IMP.Model()
universe=IMP.atom.Hierarchy.setup_particle(IMP.Particle(m, "the universe"))
# Read in the CHARMM heavy atom topology and parameter files
ff = IMP.atom.get_heavy_atom_CHARMM_parameters()
optimized_particles = []
lsc_atoms= IMP.container.ListSingletonContainer(m)

for i_prot,(pdb_file,prot_name) in enumerate((("2ptc_E.pdb","r"),("2ptc_I.pdb","l"))):
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

    rb=IMP.atom.create_rigid_body(prot)
    if i_prot > 0:
        #first body is fixed
        rb.set_coordinates_are_optimized(True)
    optimized_particles.append(rb)

    universe.add_child(prot)

    # Add non-bonded interaction (in this case, Lennard-Jones). This needs to
    # know the radii and well depths for each atom, so add them from the forcefield
    # (they can also be assigned manually using the XYZR or LennardJones
    # decorators):
    ff.add_radii(prot)
    ff.add_well_depths(prot)

    # Get a list of all atoms in the protein, and put it in a container
    atoms = IMP.atom.get_by_type(prot, IMP.atom.ATOM_TYPE)
    lsc_atoms.add_particles(atoms)

# Add a restraint for the Lennard-Jones interaction.
# Then, a LennardJonesPairScore scores a pair of atoms with the Lennard-Jones
# potential. Finally, a PairsRestraint is used which simply applies the
# LennardJonesPairScore to each pair in the ClosePairContainer.
# Cutoff is spehere surface distance for XYZR and center-to-center for XYZ decorators
nbl = IMP.container.ClosePairContainer(lsc_atoms, 0, IMP.core.RigidClosePairsFinder(), 4.0)
sf = IMP.atom.ForceSwitch(6.0, 7.0)
ps = IMP.atom.LennardJonesPairScore(sf)
m.add_restraint(IMP.container.PairsRestraint(ps, nbl))

# Finally, evaluate the score of the whole system (without derivatives)
e = m.evaluate(False)
log.info("Before optimization: {}".format(e))

out_file="2ptc.ini.tmp.pdb" #IMP.create_temporary_file("2ptc.ini.", ".pdb")
IMP.atom.write_pdb(universe,out_file)

use_simplex = False

if use_simplex:

    o = IMP.gsl.Simplex(m)
    o.set_minimum_size(.0005)
    o.set_initial_length(0.001)

    o.optimize(100)
    e = o.get_last_score()

    log.info("After simplex: {}".format(e))

# Set up optimizer
#this one needs simplex first to work on 2ptc native
#o = IMP.core.ConjugateGradients()
#This one does not need simplex
o = IMP.gsl.ConjugateGradients()
#this one is warned against in IMP as experimental, but works w/o simplex
#o = IMP.gsl.QuasiNewton()
o.set_model(m)

done=False
while not done:
    try:
        o.optimize(10)
    except IMP.base.ModelException,msg:
        log.exception()
        log.info("randomly moving after optimizer exception...")
        for d in optimized_particles:
            IMP.core.transform(d,IMP.algebra.Transformation3D(IMP.algebra.get_random_rotation_3d(IMP.algebra.get_identity_rotation_3d(),0.1),
                                                              IMP.algebra.get_random_vector_in(IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(0,0,0),
                                                                                               IMP.algebra.Vector3D(0.1,0.1,0.1)))))
    else:
        done=True

e = o.get_last_score()

log.info("After optimization {}".format(e))

out_file="2ptc.opt.tmp.pdb" #IMP.create_temporary_file("2ptc.opt.", ".pdb")
IMP.atom.write_pdb(universe,out_file)

logging.root.removeHandler(log_h_console)