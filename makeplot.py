from triplot import getsigmalevels, makesubplot2d
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import libstempo 
import matplotlib.gridspec as gridspec


gs = gridspec.GridSpec(2,2)

ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[1,:])


def makecontour(datafile, parfile, timfile, ax):
    #data = np.load('1950.after.npz')
    data = np.load(datafile)
    fitpars = list(data['fitpars'])
    res = data['res']

    #psr = libstempo.tempopulsar(parfile='before_DDGR.par', timfile='J1950+2414_T2.tim') 
    #psr = libstempo.tempopulsar(parfile='after_DDGR.par', timfile='1950.all.tim')
    psr = libstempo.tempopulsar(parfile=parfile, timfile=timfile)
    psr.fit()
    pars = psr.pars() 
    vals = psr.vals()
    errs = psr.errs()
    fitidx = [i for i,p in enumerate(pars) if p in fitpars ]
    #print fitidx
    vals0 = vals[fitidx]
    errs0 = errs[fitidx]
    parsize = vals0.size
    #print parsize
    #print 'fit for: ', fitpars

    plist = fitpars
    MarkovChain = res[:,1:]

    MCMCSize = len(MarkovChain)

    def extract(par):
        ipar = plist.index(par)
        return MarkovChain[:,ipar]*errs0[ipar] + vals0[ipar]

    M2 = extract('M2').astype(float)
    MTOT = extract('MTOT').astype(float)
    M1 = MTOT - M2
    ecc = 0.0798
    NS_MAX = 2.1
    M_SS = lambda M2: (NS_MAX - ecc * M2) / (1 + ecc)


    makesubplot2d(ax, M2, M1) 
    ax.set_xlim(0.2, 0.4)
    ax.set_ylim(1.0, 2.0)
    ax.set_xlabel(r'M$_{\rm WD}$ (M$_{\odot}$)')
    ax.set_ylabel(r'M$_{\rm MSP}$ (M$_{\odot}$)')

    ax.axhspan(1.22, 1.31, facecolor='Gray', alpha=0.3)
    ax.text(0.35, 1.23, 'RD-AIC')
    m2 = np.arange(0.2, 0.41, 0.01)
    top = np.ones(m2.size) * 2.0
    ax.fill_between(m2, M_SS(m2), top,  facecolor='g', alpha=0.3)
    ax.text(0.35, 1.93, 'QS-Nova')
    ax.axvline(0.31, color='k', linewidth=2)
    ax.arrow(0.31, 1.8, -0.01, 0, head_width=0.05, head_length=0.005, fc='k', linewidth=2)
    ax.text(0.25, 1.8, 'CB disk')

    return 

makecontour('1950.before.npz', 'before_DDGR.par', 'J1950+2414_T2.tim', ax1)
makecontour('1950.after.npz', 'after_DDGR.par', '1950.all.tim', ax2)

ax1.set_title('existing observations')
ax1.text(0.22,1.85,r'$a$',fontsize=15)
ax2.set_title('add proposed observations')
ax2.text(0.22,1.85,r'$b$',fontsize=15)

parfile = 'after_DDGR.par'
timfile = '1950.all.tim'
psr = libstempo.tempopulsar(parfile=parfile, timfile=timfile)
psr.fit()
toas = psr.toas()
res = psr.residuals() * 1.e6
errs = psr.toaerrs 
TOASIZE = toas.size
CUTMJD = 57400
oldidx = np.arange(TOASIZE)[toas < CUTMJD]
newidx = np.arange(TOASIZE)[toas > CUTMJD]
ax3.errorbar(toas[oldidx], res[oldidx], yerr=errs[oldidx], fmt='k.', ms = 1)
ax3.errorbar(toas[newidx], res[newidx], yerr=errs[newidx], fmt='r.', ms = 1)
ax3.axvline(57450, color='k', linewidth=2)
ax3.text(57450,40,'03.03.2016',color='k')

ax3.set_ylabel('Residual (ms)')
ax3.set_xlabel('MJD (day)')
ax3.text(55600,40,r'$c$',fontsize=15)

show()
