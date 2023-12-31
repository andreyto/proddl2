### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
main_pkg_name = "PRODDL_CLIENT"
script_pref = "proddl-"
pkg_subdir = "lib"

import os,sys,tempfile
from fnmatch import fnmatch

#prevent easy_install from complaining about sys.path
#sys.path.insert(0,r"@PY_INSTALL_PREFIX@")
#sys.path.insert(0,r"@PY_BUILD_INSTALL_PREFIX@")

#setup.py test command does not create entry point scripts.
#Thus, we detect when a command line looks like "* setup.py test" and modify
#it to look like "* setup.py develop --install-dir <tmp_dir> test", also
#adding <tmp_dir> to PATH and PYTHONPATH.

if len(sys.argv) >= 2 and sys.argv[-1] == "test" and sys.argv[-2] == "setup.py":

    print "Modifying path for testing"

    test_inst_dir = os.path.abspath(os.path.join("test_run","install"))
    if not os.path.exists(test_inst_dir):
        os.makedirs(test_inst_dir)
    #test_inst_dir = tempfile.mkdtemp(prefix="install.test.",suffix=".tmp",dir=test_inst_dir)
    sys.path.insert(0,test_inst_dir)
    PATH = os.environ.get("PATH",test_inst_dir)
    os.environ["PATH"] = os.pathsep.join([test_inst_dir,PATH])
    PYTHONPATH = os.environ.get("PYTHONPATH",test_inst_dir)
    os.environ["PYTHONPATH"] = os.pathsep.join([test_inst_dir,PYTHONPATH])

    pos_test = -1
    sys.argv = sys.argv[:pos_test]+["develop","--install-dir",test_inst_dir]+sys.argv[pos_test:]

#activate 'distribute', installing it if necessary
from distribute_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

def _entry_point(script_name,pkg_path):
    return script_pref+script_name+' = '+main_pkg_name+'.'+pkg_path

class PyTest(TestCommand):
    """Integrates pytest as per http://pytest.org/latest/goodpractises.html?highlight=egg
    so that you can run 'python setup.py test'"""
    def finalize_options(self):
        TestCommand.finalize_options(self)
        test_dir = os.path.join(pkg_subdir,main_pkg_name,"test")
        self.test_args = ["--verbose",test_dir]
        self.test_suite = test_dir
        #print self.test_args
        #print self.test_suite
    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

def iter_files_tree(dir,to_base=True,patt=None):
    dir_base = os.path.dirname(dir)
    if patt is not None:
        if isinstance(patt,str):
            patt = [ patt ]
    for root,dirs,files in os.walk(dir):
        for f in files:
            if patt:
                match = False
                for p in patt:
                    if fnmatch(f,p):
                        match = True
                        break
                if not match:
                    continue
            path = os.path.join(root,f)
            if to_base:
                path = os.path.relpath(path,dir_base)
            yield path

def _list_pkg_files(dir,patt=None):
    return list(iter_files_tree(os.path.join(pkg_subdir,main_pkg_name,dir),patt=patt))


setup(
    name = main_pkg_name,
    version = "1.2",
    packages = find_packages(pkg_subdir),
    package_dir = {'':pkg_subdir},
    #argh is used to auto-generate command line argument 
    #processing in entry points
    install_requires = ['argh','argcomplete','bioblend','pytest'],
	#install_requires = ['argh','argcomplete','biopython>=1.61'],
	#this will install pytest module
    tests_require=['pytest'],
    cmdclass = {'test': PyTest},
    #http://pythonhosted.org/distribute/easy_install.html#editing-and-viewing-source-packages
    #from pkg_resources import load_entry_point
    #load_entry_point('distribute==0.6.35', 'console_scripts', 'easy_install')(argv=sys.argv[1:])
	#full URL could be vcs+proto://host/path@revision#egg=project-version
	#dependency_links=["git+git://github.com/chapmanb/bcbb.git#egg=bcbb"],
    #dependency_links=["file:"+os.path.join(bcbb_src,"gff")+"#egg=bcbio_gff"],
    #dependency_links=["git+git://github.com/biopython/biopython.git@21cd969c68ca3d4f79ed5d3499766ac4eaecabb7#egg=biopython-1.61_"],
    package_data = {
        main_pkg_name: _list_pkg_files("data") + \
                _list_pkg_files("test",patt="*.py"),
    },
    entry_points = {
        #'console_scripts' is a fixed group name - it will cause
        #creation of scripts
        'console_scripts': [
            _entry_point('client','client:main'),
            ]
        },
    # metadata for upload to PyPI
    author = "Andrey Tovchigrechko",
    author_email = "andreyto@gmail.com",
    description = "Remote API Python modules for PRODDL protein-protein docking package",
    license = "GPL",
    keywords = "protein-protein docking distributed makeflow client API",
    url = "http://bitbucket.org/andreyto/proddl2", 
)

