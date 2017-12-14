#!/usr/bin/env python
import os
import numpy as np
from subprocess import call

FSIZ = 10000

wlen = 0.55
refrs = np.arange(1.3,1.8001,0.01)
refis = np.power(10.0,np.arange(-4.0,0.0001,0.1))
xs = np.power(10.0,np.arange(-4.0,2.0001,0.01))
es = np.arange(2.0,2.01,0.1)
for e in es:
    ne = '%06d'%(int(e*100000+0.5))
    if not os.path.exists(ne):
        os.makedirs(ne)
    if not os.path.isdir(ne):
        raise IOError('No such directory >>> '+ne)
for x in xs:
    nx = '%09d'%(int(x*1000000+0.5))
    for refr in refrs:
        nr = '%06d'%(int(refr*100000+0.5))
        for refi in refis:
            ni = '%06d'%(int(refi*100000+0.5))
            for e in es:
                ne = '%06d'%(int(e*100000+0.5))
                fnam = os.path.join(ne,'tmq_e%s_x%s_refr%s_refi%s.dat'%(ne,nx,nr,ni))
                if os.path.exists(fnam) and os.path.getsize(fnam) > FSIZ:
                    continue
                call('tmq_calc -b -c 4 -w %s -e %s -x %.15e -R %s -I %.15e >%s'%(wlen,e,x,refr,refi,fnam),shell=True)
