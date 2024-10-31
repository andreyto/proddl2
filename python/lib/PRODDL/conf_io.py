### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

"""Methods for loading config files"""

from . import resources
from . import json_conf
from . import util
import os

def load_config_standard_vars(config_file=None,
                              home_dir=None,
                              data_dir=None,
                              wrapper=None):

    if home_dir is None:
        home_dir = resources.get_home_dir()
    if config_file is None:
        config_file = resources.get_config_file()
    if data_dir is None:
        data_dir = resources.get_data_dir()
    if wrapper is None:
        wrapper = resources.get_wrapper()

    if util.is_string(config_file):
        return json_conf.load_config_json(config_file,
                                          variables=dict(
                                                         home_dir=home_dir,
                                                         data_dir=data_dir,
                                                         wrapper=wrapper
                                                         )
                                          )
    return config_file

def save_config(config,config_file):
    json_conf.save_config_json(config,config_file)
