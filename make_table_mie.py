#!/usr/bin/env python
import os
import numpy as np
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

DATDIR = 'mie5'
DSTDIR = 'mie'
RMIN = 1.3
RMAX = 1.8
RSTP = 0.01
REFR_MULT = 1.0e6
REFR_NDIG = 8
IMIN = -5.0
IMAX = 0.0
ISTP = 0.1
REFI_MULT = 1.0e6
REFI_NDIG = 8
XMIN = -4.0
XMAX = 4.0
XSTP = 0.001

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-d','--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('-D','--dstdir',default=DSTDIR,help='Destination directory (%default)')
parser.add_option('-r','--rmin',default=RMIN,type='float',help='Min. refractive index -- real part (%default)')
parser.add_option('-R','--rmax',default=RMAX,type='float',help='Max. refractive index -- real part (%default)')
parser.add_option('--rstp',default=RSTP,type='float',help='Refractive index step -- real part (%default)')
parser.add_option('--refr_mult',default=REFR_MULT,type='float',help='Refractive index multiplying factor -- real part (%default)')
parser.add_option('--refr_ndig',default=REFR_NDIG,type='int',help='Refractive index #digits --real part (%default)')
parser.add_option('-i','--imin',default=IMIN,type='float',help='Min. refractive index -- imaginary part, log (%default)')
parser.add_option('-I','--imax',default=IMAX,type='float',help='Max. refractive index -- imaginary part, log (%default)')
parser.add_option('--istp',default=ISTP,type='float',help='Refractive index step -- imaginary part, log (%default)')
parser.add_option('--refi_mult',default=REFI_MULT,type='float',help='Refractive index multiplying factor -- imaginary part (%default)')
parser.add_option('--refi_ndig',default=REFI_NDIG,type='int',help='Refractive index #digits --imaginary part (%default)')
parser.add_option('-x','--xmin',default=XMIN,type='float',help='Min. size parameter -- log (%default)')
parser.add_option('-X','--xmax',default=XMAX,type='float',help='Max. size parameter -- log (%default)')
parser.add_option('--xstp',default=XSTP,type='float',help='Size parameter step -- log (%default)')
(opts,args) = parser.parse_args()

if not os.path.isdir(opts.datdir):
    raise IOError('No such directory >>> '+opts.datdir)
if not os.path.exists(opts.dstdir):
    os.makedirs(opts.dstdir)
if not os.path.isdir(opts.dstdir):
    raise IOError('No such directory >>> '+opts.dstdir)

refrs = np.arange(opts.rmin,opts.rmax+0.1*opts.rstp,opts.rstp)
refis = np.power(10.0,np.arange(opts.imin,opts.imax+0.1*opts.istp,opts.istp))
fmt_refr = '%%0%dd'%(opts.refr_ndig)
fmt_refi = '%%0%dd'%(opts.refi_ndig)
for refr in refrs:
    nr = fmt_refr%(int(refr*opts.refr_mult+0.5))
    for refi in refis:
        ni = fmt_refi%(int(refi*opts.refi_mult+0.5))
        bnam = os.path.join(opts.dstdir,'mie_refr%s_refi%s'%(nr,ni))
        fnam = os.path.join(opts.datdir,'mie_refr%s_refi%s.dat'%(nr,ni))
        print fnam
        call('read_mie -M 721 -l %s -L %s -s %s -R %.22e -I %.22e <%s'%(
              opts.xmin,opts.xmax+0.1*opts.xstp,opts.xstp,refr,refi,fnam),shell=True)
        os.rename('qext.dat',bnam+'_qext.dat')
        os.rename('qsca.dat',bnam+'_qsca.dat')
        os.rename('p1.dat',bnam+'_p1.dat')
        os.rename('p2.dat',bnam+'_p2.dat')
