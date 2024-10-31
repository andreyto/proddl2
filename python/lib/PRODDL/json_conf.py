### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

"""Methods for loading config files in JSON format with variable interpolation"""

import json
import collections

def json_quote(x):
    r"""Quote a scalar for insertion into already JSON-encoded string.
    Primary use case: '["${path}/lib"]'.replace("${path}","c:\\something")
    where in JSON, Windows path should look like "c:\\\\something".
    """
    s = json.dumps((str(x),))
    assert s[:2] == '["' and s[-2:] == '"]', "Unexpected JSON structure from encoder"
    return s[2:-2]

def load_config_json(config_file,variables={}):
    with open(config_file,'r') as f:
        text = f.read()
        for key, val in list(variables.items()):
            patt = "${"+str(key)+"}"
            text = text.replace(patt,json_quote(val))
        return json.loads(text)

def save_config_json(config,config_file):
    with open(config_file,'w') as f:
        json.dump(config,f)

