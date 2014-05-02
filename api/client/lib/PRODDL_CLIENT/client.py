### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from galaxy_client import proddl

import argh
from argh import arg

from subprocess import check_call
import logging
import os, glob, shutil, json

log = logging.getLogger(__name__)

def dock(
            options,
            receptor_pdb,
            ligand_pdb,
            models_pdb,
            n_models=10
            ):
    proddl(options=options).dock(
            receptor_pdb,
            ligand_pdb,
            models_pdb,
            n_models=n_models,
            wait=True
            )

def main():
    from argh import ArghParser
    parser = ArghParser()
    parser.add_commands([dock])
    parser.dispatch()

