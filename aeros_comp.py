#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser,IndentedHelpFormatter

XMIN = 300.0
XMAX = 1300.0
XSTP = 200.0

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog input_file_names [options]')
parser.add_option('-x','--xmin',default=XMIN,type='float',help='X min. (%default)')
parser.add_option('-X','--xmax',default=XMAX,type='float',help='X max. (%default)')
parser.add_option('--xstp',default=XSTP,type='float',help='X step (%default)')
parser.add_option('--ax1_ymin',default=None,type='float',help='Axis1 Y min. (%default)')
parser.add_option('--ax1_ymax',default=None,type='float',help='Axis1 Y max. (%default)')
parser.add_option('--ax1_ystp',default=None,type='float',help='Axis1 Y step (%default)')
parser.add_option('--ax2_ymin',default=None,type='float',help='Axis1 Y min. (%default)')
parser.add_option('--ax2_ymax',default=None,type='float',help='Axis1 Y max. (%default)')
parser.add_option('--ax2_ystp',default=None,type='float',help='Axis1 Y step (%default)')
parser.add_option('--ax3_ymin',default=None,type='float',help='Axis1 Y min. (%default)')
parser.add_option('--ax3_ymax',default=None,type='float',help='Axis1 Y max. (%default)')
parser.add_option('--ax3_ystp',default=None,type='float',help='Axis1 Y step (%default)')
parser.add_option('-s','--csca',default=False,action='store_true',help='Scattering mode (%default)')
parser.add_option('-a','--cabs',default=False,action='store_true',help='Absorption mode (%default)')
parser.add_option('--xsec',default=False,action='store_true',help='Cross-section mode (%default)')
parser.add_option('--no_grid',default=False,action='store_true',help='No grid mode (%default)')
(opts,args) = parser.parse_args()

if len(args) > 0:
    fs = args
ndat = len(fs)
if ndat < 2:
    raise ValueError('Error, ndat={}'.format(ndat))

wav1 = []
aext = []
e550 = []
omeg = []
asym = []
wav2 = []
angl = []
phas = []
cnds = []
for i in range(ndat):
    w_1,aext_1,e550_1,omeg_1,asym_1 = np.loadtxt(os.path.join(fs[i],'mie_out1.dat'),unpack=True)
    wav_1,ang_1,phas_1 = np.loadtxt(os.path.join(fs[i],'mie_out2.dat'),usecols=(0,1,4),unpack=True)
    cnd_1 = (np.fabs(wav_1-550.0)<1.0e-7) & ((np.fabs(ang_1%10.0)<1.0e-7)|((ang_1<10.0)&(np.fabs(ang_1%2.5)<1.0e-7)))
    wav1.append(w_1)
    aext.append(aext_1)
    e550.append(e550_1)
    omeg.append(omeg_1)
    asym.append(asym_1)
    wav2.append(wav_1)
    angl.append(ang_1)
    phas.append(phas_1)
    cnds.append(cnd_1)

plt.interactive(True)
fig = plt.figure(1,facecolor='w',figsize=(10,6))
fig.clear()
fig.subplots_adjust(wspace=0.30,hspace=0.25)

ax1 = plt.subplot(221)
ax1.minorticks_on()
if opts.no_grid:
    ax1.grid(False)
else:
    ax1.grid(True)
if opts.csca:
    if opts.xsec:
        for i in range(ndat):
            ax1.plot(wav1[i],omeg[i]*aext[i],'o-',label=os.path.basename(fs[i]))
        ax1.set_ylabel('Scattering Cross-Section')
    else:
        for i in range(ndat):
            ax1.plot(wav1[i],omeg[i]*e550[i],'o-',label=os.path.basename(fs[i]))
        ax1.set_ylabel('Scattering Coefficient')
else:
    if opts.xsec:
        for i in range(ndat):
            ax1.plot(wav1[i],aext[i],'o-',label=os.path.basename(fs[i]))
        ax1.set_ylabel('Extinction Cross-Section')
    else:
        for i in range(ndat):
            ax1.plot(wav1[i],e550[i],'o-',label=os.path.basename(fs[i]))
        ax1.set_ylabel('Extinction Coefficient')
ax1.set_xlabel('Wavelength (nm)')
ax1.set_xlim(opts.xmin,opts.xmax)
ax1.xaxis.set_major_locator(plt.MultipleLocator(opts.xstp))
if opts.ax1_ymin is not None:
    ax1.set_ylim(bottom=opts.ax1_ymin)
if opts.ax1_ymax is not None:
    ax1.set_ylim(top=opts.ax1_ymax)
if opts.ax1_ystp is not None:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(opts.ax1_ystp))
ax1.legend(prop={'size':16},numpoints=1,ncol=ndat,loc=8,bbox_to_anchor=(1.06,1.01),frameon=False)
ax1.xaxis.set_tick_params(pad=7)
ax1.yaxis.set_label_coords(-0.12,0.5)

ax2 = plt.subplot(222)
ax2.minorticks_on()
if opts.no_grid:
    ax2.grid(False)
else:
    ax2.grid(True)
if opts.cabs:
    if opts.xsec:
        for i in range(ndat):
            ax2.plot(wav1[i],aext[i]-omeg[i]*aext[i],'o-')
        ax2.set_ylabel('Absorption Cross-Section')
    else:
        for i in range(ndat):
            ax2.plot(wav1[i],e550[i]-omeg[i]*e550[i],'o-')
        ax2.set_ylabel('Absorption Coefficient')
else:
    for i in range(ndat):
        ax2.plot(wav1[i],omeg[i],'o-')
    ax2.set_ylabel('Single Scattering Albedo')
ax2.set_xlim(opts.xmin,opts.xmax)
ax2.xaxis.set_major_locator(plt.MultipleLocator(opts.xstp))
if opts.ax2_ymin is not None:
    ax2.set_ylim(bottom=opts.ax2_ymin)
if opts.ax2_ymax is not None:
    ax2.set_ylim(top=opts.ax2_ymax)
if opts.ax2_ystp is not None:
    ax2.yaxis.set_major_locator(plt.MultipleLocator(opts.ax2_ystp))
ax2.set_xlabel('Wavelength (nm)')
ax2.xaxis.set_tick_params(pad=7)
ax2.yaxis.set_label_coords(-0.12,0.5)

ax3 = plt.subplot(223)
ax3.minorticks_on()
if opts.no_grid:
    ax3.grid(False)
else:
    ax3.grid(True)
for i in range(ndat):
    ax3.plot(wav1[i],asym[i],'o-')
ax3.set_xlabel('Wavelength (nm)')
ax3.set_ylabel('Asymmetry Parameter')
ax3.set_xlim(opts.xmin,opts.xmax)
ax3.xaxis.set_major_locator(plt.MultipleLocator(opts.xstp))
if opts.ax3_ymin is not None:
    ax3.set_ylim(bottom=opts.ax3_ymin)
if opts.ax3_ymax is not None:
    ax3.set_ylim(top=opts.ax3_ymax)
if opts.ax3_ystp is not None:
    ax3.yaxis.set_major_locator(plt.MultipleLocator(opts.ax3_ystp))
ax3.xaxis.set_tick_params(pad=7)
ax3.yaxis.set_label_coords(-0.12,0.5)

ax4 = plt.subplot(224)
ax4.minorticks_on()
if opts.no_grid:
    ax4.grid(False)
else:
    ax4.grid(True)
for i in range(ndat):
    ax4.plot(angl[i][cnds[i]],phas[i][cnds[i]],'o-')
ax4.set_xlabel('Scattering Angle (deg)')
ax4.set_ylabel('Phase Function at 550 nm')
ax4.set_xlim(0.0,180.0)
ax4.xaxis.set_major_locator(plt.MultipleLocator(30.0))
ax4.set_yscale('log')
ax4.xaxis.set_tick_params(pad=7)
ax4.yaxis.set_label_coords(-0.12,0.5)

plt.savefig('aeros_comp.pdf')
plt.draw()
