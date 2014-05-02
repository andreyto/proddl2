### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from PRODDL import resources
import pytest
from subprocess import check_call
import os, glob
from os.path import join as pjoin


pytestmark = pytest.mark.usefixtures("goto_cleandir_test")

test_data_dir = resources.get_test_data_dir()

def test_imp_import():
    import IMP


