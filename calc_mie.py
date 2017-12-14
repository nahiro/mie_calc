#!/usr/bin/env python
import os
import numpy as np
from subprocess import call

DATDIR = 'mie5'
if not os.path.exists(DATDIR):
    os.makedirs(DATDIR)
if not os.path.isdir(DATDIR):
    raise IOError('No such directory >>> '+DATDIR)

refrs = np.arange(1.3,1.8001,0.01)
refis = np.power(10.0,np.arange(-15.0,-5.0001,1.0))
xs = np.power(10.0,np.arange(-4.0,4.0001,0.001))
for refr in refrs:
    nr = '%08d'%(int(refr*1000000+0.5))
    for refi in refis:
        ni = '%018d'%(int(refi*1.0e16+0.5))
        fnam = os.path.join(DATDIR,'mie_refr%s_refi%s.dat'%(nr,ni))
        if os.path.exists(fnam):
            os.unlink(fnam)
        for x in xs:
            call('mie_calc -b -M 721 -x %.22e -R %.22e -I %.22e >>%s'%(x,refr,refi,fnam),shell=True)
