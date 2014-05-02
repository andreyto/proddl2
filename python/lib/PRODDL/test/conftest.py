### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


import pytest
import tempfile, os, shutil
from subprocess import check_call
from PRODDL import resources,util

@pytest.fixture(scope="session")
def get_test_data_dir():
    inp_test_data = resources.get_test_data_dir()
    shutil.copytree(inp_test_data,os.path.join(os.getcwd(),"test_data"))

@pytest.fixture(scope="session")
def goto_cleandir_test():
    root_test_dir = "test_run"
    if not os.path.exists(root_test_dir):
        os.makedirs(root_test_dir)
    newpath = tempfile.mkdtemp(prefix="run.test.",suffix=".tmp",dir=root_test_dir)
    os.chdir(newpath)


