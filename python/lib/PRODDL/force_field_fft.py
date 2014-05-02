### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

import config
import conf_io
import util
import imp.force_field_fft

import argh
from argh import arg

import numpy as np

import h5py
import logging
import itertools as it

log = logging.getLogger(__name__)

def write_molforce(
        receptor_pdb,
        ligand_pdb,
        molforce_file,
        home_dir=None,
        options=None):

    opt = conf_io.load_config_standard_vars(config_file=options,home_dir=home_dir)

    out = h5py.File(molforce_file,"w")
    #ljComp in FFT needs mix=1 (geometric mean for combining sigma)
    #CHARMM actually uses arithmetic mean [1](
    #http://www.charmm.org/ubbthreads/ubbthreads.php?ubb=showflat&Number=25281
    #)
    out.attrs["mix"] = np.intc(1) 

    mol_offsets = []
    pos = []
    mass = []
    eps = []
    sigma = []
    alpha = []

    ind_at = 0
    for mol_name, recs_mol in it.groupby(imp.force_field_fft.extract_molforce((receptor_pdb,ligand_pdb)),
                           lambda rec: rec["mol_name"]):
        ind_at_start = ind_at
        for rec in recs_mol:
            pos.append(rec["xyz"])
            mass.append(rec["mass"])
            eps.append(rec["eps"])
            sigma.append(rec["rad"])
            ind_at += 1
        mol_offsets.append((ind_at_start,ind_at))
    
    out["mol_offsets"] = np.asarray(mol_offsets,np.intc)

    out["pos"] = np.asarray(pos,config.T_num)
    
    out["mass"] = np.asarray(mass,config.T_num)
    
    # Our C++ code was based on OPLS and Amber LJ parametrization.
    # IMP uses CHARMM.
    # Here is how to convert the parameters:
    #
    # OPLS = 4 eps_opls [ (sigma/r)^12 - (sigma/r)^6 ] (1)
    # CHARMM = eps_charmm [ (r_min/r)^12 - 2 (r_min/r)^6 ] (2)
    # Solving for CHARMM' = 0, we see that r_min is indeed the coordinate
    # of the minimum. Ref [1] says that CHARMM parameter files actually
    # contain r_min/2. Since r_min is the equilibrium distance between
    # two atoms, r_min/2 can be viewed as a Van-der-Waals radius (e.g. 1.7A).
    # IMP returns values ~ 2.0A.
    # We find a conversion between OPLS and CHARMM parameters that will make
    # OPLS == CHARMM for all r.
    # We set r -> r_min in (1) and r -> +Inf in (2) (so that power 12 term can be
    # ignored) and solve the resulting two equations.
    # This has a single solution where 
    # eps_opls = eps_charmm 
    # and 
    # sigma = 2^(-1/6) r_min
    # We also invert sign for os eps in our implementation (IMP returns negative values)
    # eps_opls returned by IMP is still about twice as large.
    # CHARMM uses these units (AKMA: Angstroms, Kilocalories/Mole, Atomic mass units)
    # http://www.charmm.org/documentation/c34b1/usage.html
    # Our C++ code was getting parameters from MMTK, which used internal units of kJ/mol
    # https://bitbucket.org/khinsen/mmtk/src/a5fb7cde6f5b978cafceab1dcad9e1d7984f4e90/MMTK/Units.py
    # So, we need to multiply eps_charmm by 4.184.
    
    out["eps"] = - 4.184 * np.asarray(eps,config.T_num)
    
    # CHARMM parameter is 1/2 * r_min [1]
    out["sigma"] = 2*2**(-1./6)*np.asarray(sigma,config.T_num)
    
    alpha = np.zeros(len(sigma),config.T_num)
    alpha[:] = opt["scan"]["alpha"]
    out["alpha"] = alpha

    out.close()
