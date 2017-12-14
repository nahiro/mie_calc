#!/usr/bin/env python
import os
import sys
import struct
import numpy as np
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

fmt_data = '=ddddddiii'
dsiz = struct.calcsize(fmt_data)

# Default values
DATDIR = '.'
DSTDIR = 'tmq'
NDGS = 4
NANG = 361
WLEN = 0.55
RMIN = 1.3
RMAX = 1.8
RSTP = 0.01
REFR_MULT = 1.0e5
REFR_NDIG = 6
IMIN = -4.0
IMAX = 0.0
ISTP = 0.1
REFI_MULT = 1.0e5
REFI_NDIG = 6
XMIN = -4.0
XMAX = 2.0
XSTP = 0.01
X_MULT = 1.0e6
X_NDIG = 9
EMIN = 2.0
EMAX = 2.0
ESTP = 0.1
E_MULT = 1.0e5
E_NDIG = 6
ELST = []

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-d','--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('-D','--dstdir',default=DSTDIR,help='Destination directory (%default)')
parser.add_option('-c','--ndgs',default=NDGS,type='int',help='Control# (%default)')
parser.add_option('-M','--nang',default=NANG,type='int',help='#Angles (%default)')
parser.add_option('-w','--wlen',default=WLEN,type='float',help='Wavelength in um (%default)')
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
parser.add_option('--x_mult',default=X_MULT,type='float',help='Size parameter multiplying factor (%default)')
parser.add_option('--x_ndig',default=X_NDIG,type='int',help='Size parameter #digits (%default)')
parser.add_option('-e','--emin',default=EMIN,type='float',help='Min. aspect ratio (%default)')
parser.add_option('-E','--emax',default=EMAX,type='float',help='Max. aspect ratio (%default)')
parser.add_option('--estp',default=ESTP,type='float',help='Aspect ratio step (%default)')
parser.add_option('--reps',default=None,type='float',help='Actual aspect ratio (%default)')
parser.add_option('--e_mult',default=E_MULT,type='float',help='Aspect ratio multiplying factor (%default)')
parser.add_option('--e_ndig',default=E_NDIG,type='int',help='Aspect ratio #digits (%default)')
parser.add_option('-l','--elst',default=ELST,action='append',help='Exception list (%default)')
parser.add_option('-S','--skip_error',default=False,action='store_true',help='Skip errors (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
(opts,args) = parser.parse_args()

if not os.path.exists(opts.dstdir):
    os.makedirs(opts.dstdir)
if not os.path.isdir(opts.dstdir):
    raise IOError('No such directory >>> '+opts.dstdir)

refrs = np.arange(opts.rmin,opts.rmax+0.1*opts.rstp,opts.rstp)
refis = np.power(10.0,np.arange(opts.imin,opts.imax+0.1*opts.istp,opts.istp))
xs = np.power(10.0,np.arange(opts.xmin,opts.xmax+0.1*opts.xstp,opts.xstp))
es = np.arange(opts.emin,opts.emax+0.1*opts.estp,opts.estp)
fmt_refr = '{{:0{:d}d}}'.format(opts.refr_ndig)
fmt_refi = '{{:0{:d}d}}'.format(opts.refi_ndig)
fmt_x = '{{:0{:d}d}}'.format(opts.x_ndig)
fmt_e = '{{:0{:d}d}}'.format(opts.e_ndig)
for e in es:
    ne = fmt_e.format(int(e*opts.e_mult+0.5))
    ne2 = '{:02d}'.format(int(e*10.0+0.5))
    if opts.reps is not None:
        reps = opts.reps
        na = fmt_e.format(int(opts.reps*opts.e_mult+0.5))
    else:
        reps = e
        na = ne
    for refr in refrs:
        nr = fmt_refr.format(int(refr*opts.refr_mult+0.5))
        for refi in refis:
            ni = fmt_refi.format(int(refi*opts.refi_mult+0.5))
            bnam = os.path.join(opts.dstdir,'tmq_e{}_refr{}_refi{}'.format(na,nr,ni))
            gnam = 'dummy.dat'
            flag = False
            xmin = xs[0]
            xmax = None
            if os.path.exists(bnam+'.dat'):
                print bnam,'exists'
                continue
            print bnam
            with open('temp.dat','w') as fp:
                for indx,x in enumerate(xs):
                    nx = fmt_x.format(int(x*opts.x_mult+0.5))
                    fnam = os.path.join(opts.datdir,ne2,ne,'tmq_e{}_x{}_refr{}_refi{}.dat'.format(ne,nx,nr,ni))
                    try:
                        with open(fnam,'r') as gp:
                            data = struct.unpack(fmt_data,gp.read(dsiz))
                        ierr = data[8]
                        if ierr != 0:
                            if opts.skip_error or fnam in opts.elst:
                                raise ValueError('Error, ierr={}'.format(ierr))
                            else:
                                raise UserWarning('Error, ierr={}'.format(ierr))
                        if flag:
                            if opts.skip_error or fnam in opts.elst:
                                raise ValueError('Flag error >>> '+fnam)
                            else:
                                raise UserWarning('Flag error >>> '+fnam)
                    except (IOError, ValueError, struct.error) as inst:
                        if opts.verbose:
                            print inst,' >>> '+fnam
                        if not flag:
                            xmax = xs[indx-1]
                        flag = True
                        fnam = gnam
                    if not os.path.exists(fnam):
                        raise IOError('No such file >>> {}\n'.format(fnam))
                    gnam = fnam
                    with open(fnam,'r') as gp:
                        data = gp.read()
                    fp.write(data)
                if xmax is None:
                    xmax = xs[-1]
            retcode = call('read_tmq -c {} -w {} -M {} -l {} -L {} -s {} -e {} -R {} -I {:.15e} <temp.dat'.format(
                            opts.ndgs,opts.wlen,opts.nang,opts.xmin,opts.xmax+0.1*opts.xstp,opts.xstp,reps,refr,refi),shell=True)
            if retcode != 0:
                raise ValueError('Error in read_tmq.')
            os.rename('qext.dat',bnam+'_qext.dat')
            os.rename('qsca.dat',bnam+'_qsca.dat')
            os.rename('p1.dat',bnam+'_p1.dat')
            os.rename('p2.dat',bnam+'_p2.dat')
            with open(bnam+'.dat','w') as fp:
                fp.write('{:10.6f} {:10.6f} {:10.6f}\n'.format(np.log10(xmin),np.log10(xmax),opts.xstp))
os.remove('temp.dat')
