#!/usr/bin/env python
import sys
import numpy as np
import struct

if len(sys.argv) < 2:
    sys.stderr.write('Usage: read_tmq.py file_name\n')
    sys.exit()
with open(sys.argv[1],'r') as fp:
    fmt = '='
    fmt += 'd' # wlen
    fmt += 'd' # xval
    fmt += 'd' # refr
    fmt += 'd' # refi
    fmt += 'd' # tmx_reps
    fmt += 'd' # tmx_delt
    fmt += 'i' # tmx_ndgs
    fmt += 'i' # tmx_shap
    fmt += 'i' # ierr
    fmt += 'i' # dummy
    fmt += 'd' # reff
    fmt += 'd' # veff
    fmt += 'd' # cext
    fmt += 'd' # csca
    fmt += 'd' # walb
    fmt += 'd' # gsym
    fmt += 'i' # lmax
    data = struct.unpack(fmt,fp.read(struct.calcsize(fmt)))
    wlen = data[0]
    xval = data[1]
    refr = data[2]
    refi = data[3]
    tmx_reps = data[4]
    tmx_delt = data[5]
    tmx_ndgs = data[6]
    tmx_shap = data[7]
    ierr = data[8]
    dummy = data[9]
    reff = data[10]
    veff = data[11]
    cext = data[12]
    csca = data[13]
    walb = data[14]
    gsym = data[15]
    lmax = data[16]
    fmt = '=%dd'%(lmax*6)
    data = np.array(struct.unpack(fmt,fp.read(struct.calcsize(fmt)))).reshape(6,-1)
    alp1 = data[0]
    alp2 = data[1]
    alp3 = data[2]
    alp4 = data[3]
    bet1 = data[4]
    bet2 = data[5]
    fmt = '=i'
    data = struct.unpack(fmt,fp.read(struct.calcsize(fmt)))
    tmx_nang = data[0]
    fmt = '=%dd'%(tmx_nang*6)
    data = np.array(struct.unpack(fmt,fp.read(struct.calcsize(fmt)))).reshape(6,-1)
    f11 = data[0]
    f22 = data[1]
    f33 = data[2]
    f44 = data[3]
    f12 = data[4]
    f34 = data[5]

sys.stdout.write('WLEN: %13.6e\n'%(wlen))
sys.stdout.write('XVAL: %13.6e\n'%(xval))
sys.stdout.write('REFR: %13.6e\n'%(refr))
sys.stdout.write('REFI: %13.6e\n'%(refi))
sys.stdout.write('REPS: %13.6e\n'%(tmx_reps))
sys.stdout.write('DELT: %13.6e\n'%(tmx_delt))
sys.stdout.write('NDGS: %13d\n'%(tmx_ndgs))
sys.stdout.write('SHAP: %13d\n'%(tmx_shap))
sys.stdout.write('IERR: %13d\n'%(ierr))
sys.stdout.write('REFF: %13.6e\n'%(reff))
sys.stdout.write('VEFF: %13.6e\n'%(veff))
sys.stdout.write('CEXT: %13.6e\n'%(cext))
sys.stdout.write('CSCA: %13.6e\n'%(csca))
sys.stdout.write('WALB: %13.6e\n'%(walb))
sys.stdout.write('GSYM: %13.6e\n'%(gsym))
sys.stdout.write('LMAX: %13d\n'%(lmax))
imax = lmax-1
sys.stdout.write('ALP1:\n')
for i in range(lmax):
    sys.stdout.write('%13.6e%s'%(alp1[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('ALP2:\n')
for i in range(lmax):
    sys.stdout.write('%13.6e%s'%(alp2[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('ALP3:\n')
for i in range(lmax):
    sys.stdout.write('%13.6e%s'%(alp3[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('ALP4:\n')
for i in range(lmax):
    sys.stdout.write('%13.6e%s'%(alp4[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('BET1:\n')
for i in range(lmax):
    sys.stdout.write('%13.6e%s'%(bet1[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('BET2:\n')
for i in range(lmax):
    sys.stdout.write('%13.6e%s'%(bet2[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('NANG: %13d\n'%(tmx_nang))
imax = tmx_nang-1
sys.stdout.write('F11:\n')
for i in range(tmx_nang):
    sys.stdout.write('%13.6e%s'%(f11[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('F22:\n')
for i in range(tmx_nang):
    sys.stdout.write('%13.6e%s'%(f22[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('F33:\n')
for i in range(tmx_nang):
    sys.stdout.write('%13.6e%s'%(f33[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('F44:\n')
for i in range(tmx_nang):
    sys.stdout.write('%13.6e%s'%(f44[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('F12:\n')
for i in range(tmx_nang):
    sys.stdout.write('%13.6e%s'%(f12[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('F34:\n')
for i in range(tmx_nang):
    sys.stdout.write('%13.6e%s'%(f34[i],'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
