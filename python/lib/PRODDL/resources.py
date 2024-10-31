### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Configuration and resources of this package"""
from . import config
from pkg_resources import Requirement, resource_filename, resource_string
import os, sys

pkg_name = config.pkg_name
pkg_data = pkg_name+"/data"

def get_pkg_data_dir():
    return resource_filename(Requirement.parse(pkg_name),pkg_data)

def get_pkg_data_file(name):
    return resource_filename(Requirement.parse(pkg_name),pkg_data+"/"+name)

def get_pkg_data_string(name):
    return resource_string(Requirement.parse(pkg_name),pkg_data+"/"+name)

def get_data_file(name=None,name_pkg=None):
    if name is None:
        return get_pkg_data_file(name_pkg)
    else:
        return name

def get_data_string(name=None,name_pkg=None):
    if name is None:
        return get_pkg_data_string(name_pkg)
    else:
        with open(name,'r') as f:
            return f.read()

def set_data_string(s,name):
    with open(name,'w') as f:
        f.write(s)

def get_default_conf_file(ext=".json"):
    conf_file = os.environ.get(pkg_name+"_CONF")
    if not conf_file:
        #assume this is called from a script and put conf file in the same dir
        conf_file = os.path.join(os.path.dirname(sys.argv[0]),pkg_name+ext)
    return conf_file

def get_home_dir():
    return os.environ.get(config.prefix_env+"_HOME",None)

def get_config_file():
    return os.environ.get(config.prefix_env+"_CONFIG",None)

def get_data_dir():
    return os.environ.get(config.prefix_env+"_DATA_DIR",None)

def get_test_data_dir():
    return os.environ.get(config.prefix_env+"_TEST_DATA_DIR",None)

def get_wrapper():
    return os.environ.get(config.prefix_env+"_WRAPPER",None)
