### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

import json
from collections import OrderedDict

def gen_options():
        dockParams = OrderedDict()
        dockParams.setdefault('recId','r')
        dockParams.setdefault('ligId','l')        
        dockParams.setdefault('chainIdR','*')
        dockParams.setdefault('chainIdL','*')
        dockParams.setdefault('alpha',0.4)
        dockParams.setdefault('zapFlexResidues',False)
        dockParams.setdefault("gridStep",0.15)
        dockParams.setdefault("cutOffFft",0.9)
        dockParams.setdefault("angleStepDeg",10)
        dockParams.setdefault("anglesDir","@PRODDL_ANGLES_DIR@")
        dockParams.setdefault('alphaFlexResidues',dockParams["alpha"])
        dockParams.setdefault('maxNRigidMatches',30000)
        dockParams.setdefault('maxNTrans',10)
        dockParams.setdefault('maxNTransInp',10000)
        dockParams.setdefault('maxValCorr',-10.)
        dockParams.setdefault('doClusterTranslations',False)
        dockParams.setdefault('clusterRadiusTrans',0.5)
        dockParams.setdefault('nMultimer',0)
        dockParams.setdefault('maxRmsdSymm',0.8)
        dockParams.setdefault('logLevel',1)
        dockParams.setdefault('testMode',0)
        dockParams.setdefault('testMaxRot',1)
        dockParams.setdefault('mmForceField','charmm')
        dockParams.setdefault('potentialName','ljComp')        
        dockParams.setdefault('lj',{})
        lj = dockParams["lj"]
        lj.setdefault('uniform',{})
        uniform = lj["uniform"]
        uniform.setdefault('sigma',0.33)
        uniform.setdefault('eps',0.46) #0.3 #0.46
        uniform.setdefault('method','scaled')
        return dockParams

fft_opt = gen_options()
opt = dict(scan=fft_opt)
print json.dumps(opt,indent=4)
