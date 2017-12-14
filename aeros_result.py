#!/usr/bin/env python
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import splrep,splev,interp1d
from StringIO import StringIO
from optparse import OptionParser,IndentedHelpFormatter

# Constants
DSR_AOD_ID = 0
DSR_AER_ID = 1
DSR_H2O_ID = 2
AUR_ID = 9
SSR_ID = 13

# Default values
FNAM = 'nohup.out'
FIGNAM = 'aeros_result.pdf'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-F','--fnam',default=FNAM,help='File name (%default)')
parser.add_option('-O','--fignam',default=FIGNAM,help='Figure name (%default)')
parser.add_option('-t','--title',default=None,help='Figure title (%default)')
parser.add_option('--dsr_min',default=None,type='float',help='DSR min in mW/m2/nm (%default)')
parser.add_option('--dsr_max',default=None,type='float',help='DSR max in mW/m2/nm (%default)')
parser.add_option('--aur_min',default=None,type='float',help='AUR min in mW/m2/nm/sr (%default)')
parser.add_option('--aur_max',default=None,type='float',help='AUR max in mW/m2/nm/sr (%default)')
parser.add_option('--ssr_min',default=None,type='float',help='SSR min in mW/m2/nm/sr (%default)')
parser.add_option('--ssr_max',default=None,type='float',help='SSR max in mW/m2/nm/sr (%default)')
parser.add_option('-S','--skip',default=False,action='store_true',help='Skip mode (%default)')
parser.add_option('-b','--batch',default=False,action='store_true',help='Batch mode (%default)')
(opts,args) = parser.parse_args()

lines = None
with open(opts.fnam,'r') as fp:
    lines = fp.readlines()
nline = len(lines)

i = 0
for i in range(nline):
    if re.search('3rd MINUIT results',lines[i]):
        break
indx = i
n = 0
for i in range(indx,nline):
    item = lines[i].split()
    nitem = len(item)
    if n==3 and nitem==8:
        break
    n = nitem
indx = i

text =''
for i in range(indx,nline):
    item = lines[i].split()
    nitem = len(item)
    if nitem != 8:
        break
    text += lines[i]

rid,pid,th_los,ph_los,wlen,data,dsim,derr = np.loadtxt(StringIO(text),unpack=True)

rids = np.array(sorted(list(set(rid))))
pids = np.array(sorted(list(set(pid))))

if not opts.batch:
    plt.interactive(True)
fig = plt.figure(1)
fig.set_facecolor('w')
fig.set_size_inches((8.0,6.0),forward=True)
plt.subplots_adjust(top=0.87)
pdf = PdfPages(opts.fignam)
for i in rids:
    for j in pids:
        cnd = (rid==i) & (pid==j)
        if not np.any(cnd):
            continue
        if cnd.sum() < 2:
            continue
        wlen_cnd = wlen[cnd]
        data_cnd = data[cnd]
        dsim_cnd = dsim[cnd]
        derr_cnd = derr[cnd]
        th_cnd = th_los[cnd]
        ph_cnd = ph_los[cnd]
        fig.clear()
        ax1 = plt.subplot(111)
        ax1.minorticks_on()
        ax1.grid(True)
        if i == DSR_H2O_ID:
            indx = np.argmin(wlen_cnd)
            x1 = wlen_cnd[indx]
            y1 = data_cnd[indx]/dsim_cnd[indx]
            indx = np.argmax(wlen_cnd)
            x2 = wlen_cnd[indx]
            y2 = data_cnd[indx]/dsim_cnd[indx]
            f = interp1d([x1,x2],[y1,y2])
            dsim_cnd *= f(wlen_cnd)
        else:
            indx = np.argmin(wlen_cnd)
            x1 = wlen_cnd[indx]
            indx = np.argmax(wlen_cnd)
            x2 = wlen_cnd[indx]
        x = np.linspace(x1,x2,100)
        y1 = splev(x,splrep(wlen_cnd,data_cnd))
        y2 = splev(x,splrep(wlen_cnd,dsim_cnd))
        ax1.plot(x,y1,'b-')
        ax1.plot(wlen_cnd,data_cnd,'bo',label='Data')
        ax1.errorbar(wlen_cnd,data_cnd,yerr=derr_cnd,color='b',fmt=None)
        ax1.plot(x,y2,'r:')
        ax1.plot(wlen_cnd,dsim_cnd,'rd',label='Simulation')
        ax1.set_xlabel('Wavelength (nm)')
        if i == SSR_ID:
            ax1.set_ylabel('Radiance (mW/m$^{2}$/nm/sr)')
            ax1.set_title(r'SSR ($\theta$ = %.2f, $\phi$ = %.2f)'%(th_cnd.mean(),ph_cnd.mean()))
            if opts.ssr_min is not None:
                ax1.set_ylim(bottom=opts.ssr_min)
            if opts.ssr_max is not None:
                ax1.set_ylim(top=opts.ssr_max)
        elif i == AUR_ID:
            ax1.set_ylabel('Radiance (mW/m$^{2}$/nm/sr)')
            ax1.set_title(r'AUR ($\theta$ = %.2f, $\phi$ = %.2f)'%(th_cnd.mean(),ph_cnd.mean()))
            if opts.aur_min is not None:
                ax1.set_ylim(bottom=opts.aur_min)
            if opts.aur_max is not None:
                ax1.set_ylim(top=opts.aur_max)
        else:
            ax1.set_ylabel('Irradiance (mW/m$^{2}$/nm)')
            ax1.set_title(r'DSR ($\theta$ = %.2f, $\phi$ = %.2f)'%(th_cnd.mean(),ph_cnd.mean()))
            if opts.dsr_min is not None:
                ax1.set_ylim(bottom=opts.dsr_min)
            if opts.dsr_max is not None:
                ax1.set_ylim(top=opts.dsr_max)
        ax1.xaxis.set_tick_params(pad=7)
        ax1.yaxis.set_label_coords(-0.10,0.5)
        ax1.legend()
        if opts.title is not None:
            plt.suptitle(opts.title,size=16,y=0.97)
        plt.savefig(pdf,format='pdf')
        if not opts.batch:
            plt.draw()
            if not opts.skip:
                sys.stderr.write('Type \'q\' to exit.')
                line = sys.stdin.readline()
                if re.search('q',line):
                    pdf.close()
                    sys.exit(0)
pdf.close()
