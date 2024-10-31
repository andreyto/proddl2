### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from . import resources
from . import conf_io
from . import util
from . import force_field_fft

from MGT.BatchMakeflow import makeflow

import argh
from argh import arg

from subprocess import check_call
import logging
import os, glob, shutil, tempfile

log = logging.getLogger(__name__)

def dock(
        receptor_pdb,
        ligand_pdb,
        model_pdb,
        n_models=10,
        home_dir=None,
        options=None,
        makeflow_args="",
        web=False,
        run_no=False,
        test_mode=False
        ):
    
    opt = conf_io.load_config_standard_vars(config_file=options,home_dir=home_dir)
    
    workflow = "dock_top.mkf"
    molforce_file = "molforce.dat"
    scan_opt_file = "scan_opt.json"
    rot_scan_list = "rot_scan_list.tab"
    res_file = "res.dat"

    wrapper = opt["wrapper"]

    with makeflow(
        makeflow_bin=opt["makeflow_bin"],
        wrapper=wrapper,
        workflow=workflow,
        makeflow_args=makeflow_args,
        workflow_script=None,
        web=web,
        run=not run_no
        ) as mf_top:

        cmd = """\
        {wrapper} \
        proddl-dock write-molforce \
        {receptor_pdb} \
        {ligand_pdb} \
        {molforce_file}""".format(**locals())

        mf_top.task(
                cmd=cmd,
                targets=[molforce_file],
                inputs=[receptor_pdb,ligand_pdb,options],
                is_local=False
                )

        opt_scan = opt["scan"]
        
        if test_mode:
            opt_scan["testMode"] = 1

        conf_io.save_config(opt_scan,scan_opt_file)

        n_ang = -1 #first line is a header
        with open(opt_scan["anglesFile"],"r") as ang:
            for l in ang:
                n_ang += 1

        if test_mode:
            n_ang = min(opt_scan["testMaxRot"],n_ang)

        n_scans = min(1000,n_ang)
        n_ang_scan = int(n_ang/n_scans)

        scan_res_files = []
        
        for start_scan in range(0,n_ang,n_ang_scan):
            
            end_scan = min(start_scan+n_ang_scan,n_ang)
            scan_res_file = "scan_res.{:04}-{:04}.dat".format(start_scan,end_scan)

            scan_res_files.append(scan_res_file)
            
            cmd = """\
            {wrapper} \
            proddl-dock-fft \
            --options {scan_opt_file} \
            --task rot-scan \
            --fft-rot-grid-start {start_scan} \
            --fft-rot-grid-end {end_scan} \
            --molforce-params {molforce_file} \
            --fft-rot-scan-res {scan_res_file}
            """.format(**locals())
            
            mf_top.task(
                    cmd=cmd,
                    targets=[scan_res_file],
                    inputs=[molforce_file,scan_opt_file]
                    )

        with open(rot_scan_list,"w") as out:
            out.write("\n".join(scan_res_files)+"\n")

        cmd = """\
        {wrapper} \
        proddl-dock-fft \
        --options {scan_opt_file} \
        --task gather \
        --fft-rot-scan-list {rot_scan_list} \
        --fft-res {res_file} \
        --molforce-params {molforce_file}
        """.format(**locals())
        
        mf_top.task(
                cmd=cmd,
                targets=[res_file],
                inputs=[molforce_file,scan_opt_file,rot_scan_list]+scan_res_files,
                is_local=False
                )

        cmd = """\
        {wrapper} \
        proddl-export \
        --model-inp {res_file} \
        --model-ind-start 0 \
        --model-ind-end {n_models} \
        --pdb-inp-rec {receptor_pdb} \
        --pdb-inp-lig {ligand_pdb} \
        --format-out pdb_nmr \
        --model-out {model_pdb}
        """.format(**locals())
        
        mf_top.task(
                cmd=cmd,
                targets=[model_pdb],
                inputs=[res_file,receptor_pdb,ligand_pdb],
                is_local=False
                )

def main():
    from argh import ArghParser
    parser = ArghParser()
    parser.add_commands([dock,force_field_fft.write_molforce])
    parser.dispatch()

