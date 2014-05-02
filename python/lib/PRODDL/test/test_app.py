### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from PRODDL import resources, conf_io
import pytest
from subprocess import check_call
import os, glob
from os.path import join as pjoin


pytestmark = pytest.mark.usefixtures("goto_cleandir_test")

test_data_dir = resources.get_test_data_dir()
opt = conf_io.load_config_standard_vars()

pdb_inp_rec = "{test_data_dir}/pdb/2ptc_E.pdb".format(test_data_dir=test_data_dir)
pdb_inp_lig = "{test_data_dir}/pdb/2ptc_I.pdb".format(test_data_dir=test_data_dir)

molforce_file = "molforce.dat"
scan_opt_file = "scan_opt.json"
scan_res_file = "scan_res.dat"
rot_scan_list = "rot_scan_list.tab"
res_file = "res.dat"
export_pdb_file = "res.pdb"

wrapper = opt["wrapper"]

def no_test_molforce():
    cmd = """\
    {wrapper} \
    proddl-dock write-molforce \
	{pdb_inp_rec} \
	{pdb_inp_lig} \
	{molforce_file}""".format(
                       wrapper=opt["wrapper"],
                       test_data_dir=test_data_dir,
                       molforce_file=molforce_file,
                       pdb_inp_rec=pdb_inp_rec,
                       pdb_inp_lig=pdb_inp_lig)
    print "opt=",opt
    print "Executing command: {}".format(cmd)
    check_call(cmd,shell=True)

def no_test_c_fft_scan_slave():
    no_test_molforce()
    conf_io.save_config(opt["scan"],scan_opt_file)
    cmd = """\
    {wrapper} \
    proddl-dock-fft \
    --options {scan_opt_file} \
    --task rot-scan \
    --fft-rot-grid-start 0 \
    --fft-rot-grid-end 2 \
    --molforce-params {molforce_file} \
    --fft-rot-scan-res {scan_res_file}
    """.format(
               wrapper=opt["wrapper"],
               molforce_file=molforce_file,
               scan_opt_file=scan_opt_file,
               scan_res_file=scan_res_file
               )
    print "Executing command: {}".format(cmd)
    check_call(cmd,shell=True)

def no_test_c_fft_scan_master():
    no_test_c_fft_scan_slave()
    with open(rot_scan_list,"w") as out:
        out.write(scan_res_file+"\n")
    cmd = """\
    {wrapper} \
    proddl-dock-fft \
    --options {scan_opt_file} \
    --task gather \
    --fft-rot-scan-list {rot_scan_list} \
    --fft-res {res_file} \
    --molforce-params {molforce_file}
    """.format(
               wrapper=opt["wrapper"],
               molforce_file=molforce_file,
               scan_opt_file=scan_opt_file,
               rot_scan_list=rot_scan_list,
               res_file=res_file
               )
    print "Executing command: {}".format(cmd)
    check_call(cmd,shell=True)

def no_test_c_export():
    no_test_c_fft_scan_master()
    cmd = """\
    {wrapper} \
    proddl-export \
    --model-inp {res_file} \
    --model-ind-start 0 \
    --model-ind-end 10 \
    --pdb-inp-rec {pdb_inp_rec} \
    --pdb-inp-lig {pdb_inp_lig} \
    --format-out pdb_nmr \
    --model-out {export_pdb_file}
    """.format(
               wrapper=opt["wrapper"],
               molforce_file=molforce_file,
               rot_scan_list=rot_scan_list,
               res_file=res_file,
               pdb_inp_rec=pdb_inp_rec,
               pdb_inp_lig=pdb_inp_lig,
               export_pdb_file=export_pdb_file
               )
    print "Executing command: {}".format(cmd)
    check_call(cmd,shell=True)

def test_dock():
    cmd = """\
    {wrapper} \
    proddl-dock \
    dock \
    {pdb_inp_rec} \
    {pdb_inp_lig} \
    {export_pdb_file} \
    --makeflow-args "-T local" \
    --test-mode 
    """.format(**globals())
    print "Executing command: {}".format(cmd)
    check_call(cmd,shell=True)

