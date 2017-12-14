#!/usr/bin/env python
import sys
import numpy as np
import struct

if len(sys.argv) < 2:
    sys.stderr.write('Usage: read_mie.py file_name\n')
    sys.exit()
with open(sys.argv[1],'r') as fp:
    fmt = '='
    fmt += 'd' # mie_xval
    fmt += 'd' # mie_refr
    fmt += 'd' # mie_refi
    fmt += 'd' # mie_icut
    fmt += 'i' # PERFCT
    fmt += 'i' # ANYANG
    fmt += 'd' # QEXT
    fmt += 'd' # QSCA
    fmt += 'd' # GQSC
    fmt += '2d' # SFORW
    fmt += '2d' # SBACK
    fmt += '4d' # TFORW
    fmt += '4d' # TBACK
    fmt += 'd' # SPIKE
    fmt += 'i' # mie_nmom
    fmt += 'i' # mie_dmom
    fmt += 'i' # mie_ipol
    data = struct.unpack(fmt,fp.read(struct.calcsize(fmt)))
    mie_xval = data[0]
    mie_refr = data[1]
    mie_refi = data[2]
    mie_icut = data[3]
    PERFCT = data[4]
    ANYANG = data[5]
    QEXT = data[6]
    QSCA = data[7]
    GQSC = data[8]
    SFORW = data[9]+data[10]*1.0j
    SBACK = data[11]+data[12]*1.0j
    TFORW = np.array([data[13]+data[14]*1.0j,data[15]+data[16]*1.0j])
    TBACK = np.array([data[17]+data[18]*1.0j,data[19]+data[20]*1.0j])
    SPIKE = data[21]
    mie_nmom = data[22]
    mie_dmom = data[23]
    mie_ipol = data[24]
    if mie_nmom > 0:
        fmt = '=%dd'%((mie_nmom+1)*4)
        PMOM = np.array(struct.unpack(fmt,fp.read(struct.calcsize(fmt)))).reshape(4,-1)
    fmt = '=i'
    data = struct.unpack(fmt,fp.read(struct.calcsize(fmt)))
    mie_nang = data[0]
    fmt = '=%dd'%(mie_nang*2)
    data = np.array(struct.unpack(fmt,fp.read(struct.calcsize(fmt)))).reshape(-1,2)
    S1 = data[:,0]+data[:,1]*1.0j
    data = np.array(struct.unpack(fmt,fp.read(struct.calcsize(fmt)))).reshape(-1,2)
    S2 = data[:,0]+data[:,1]*1.0j

sys.stdout.write('XVAL: %22.15e\n'%(mie_xval))
sys.stdout.write('REFR: %22.15e\n'%(mie_refr))
sys.stdout.write('REFI: %22.15e\n'%(mie_refi))
sys.stdout.write('ICUT: %22.15e\n'%(mie_icut))
sys.stdout.write('PFCT: %22d\n'%(PERFCT))
sys.stdout.write('ANYA: %22d\n'%(ANYANG))
sys.stdout.write('QEXT: %22.15e\n'%(QEXT))
sys.stdout.write('QSCA: %22.15e\n'%(QSCA))
sys.stdout.write('GQSC: %22.15e\n'%(GQSC))
sys.stdout.write('SFRW: %22.15e %22.15e\n'%(SFORW.real,SFORW.imag))
sys.stdout.write('SBCK: %22.15e %22.15e\n'%(SBACK.real,SBACK.imag))
sys.stdout.write('TFW1: %22.15e %22.15e\n'%(TFORW[0].real,TFORW[0].imag))
sys.stdout.write('TFW2: %22.15e %22.15e\n'%(TFORW[1].real,TFORW[1].imag))
sys.stdout.write('TBK1: %22.15e %22.15e\n'%(TBACK[0].real,TBACK[0].imag))
sys.stdout.write('TBK2: %22.15e %22.15e\n'%(TBACK[1].real,TBACK[1].imag))
sys.stdout.write('SPIK: %22.15e\n'%(SPIKE))
sys.stdout.write('NMOM: %22d\n'%(mie_nmom))
sys.stdout.write('DMOM: %22d\n'%(mie_dmom))
sys.stdout.write('IPOL: %22d\n'%(mie_ipol))
if mie_nmom > 0:
    for j in range(4):
        sys.stdout.write('PMOM%d:\n'%(j+1))
        for i in range(mie_nmom+1):
            sys.stdout.write('%22.15e%s'%(PMOM[j][i],'\n' if i%8 == 7 else ('\n' if i == mie_nmom else ' ')))
sys.stdout.write('NANG: %22d\n'%(mie_nang))
imax = mie_nang-1
sys.stdout.write('S1(real):\n')
for i in range(mie_nang):
    sys.stdout.write('%22.15e%s'%(S1[i].real,'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('S1(imag):\n')
for i in range(mie_nang):
    sys.stdout.write('%22.15e%s'%(S1[i].imag,'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('S2(real):\n')
for i in range(mie_nang):
    sys.stdout.write('%22.15e%s'%(S2[i].real,'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
sys.stdout.write('S2(imag):\n')
for i in range(mie_nang):
    sys.stdout.write('%22.15e%s'%(S2[i].imag,'\n' if i%8 == 7 else ('\n' if i == imax else ' ')))
