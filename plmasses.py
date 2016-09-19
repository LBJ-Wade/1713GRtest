from triplot import getsigmalevels, makesubplot2d
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import libstempo 

#def M1(pf):
    #Tsun = Decimal('4.925490947')*Decimal('0.000001') #Tsun == GM/c^3 in seconds
    #m2 = pf.M2[0]
    #Pb = pf.PB[0]*secperday
    #a = pf.A1[0]
    #if pf.__dict__.has_key('KIN'):
        #I = pf.KIN[0]/180*PI
        #sini = Decimal(sin(float(I)))
    #else:
        #sini = pf.SINI[0]
    ##result = sqrt(930.998*m2**3*Pb**2/a**3) - m2
    ##return (Pb/2/PI*Decimal(sqrt(G*(m2*Decimal(str(sini)))**3/a**3))-m2)/Msun
    #return Pb/2/PI*((Tsun*(m2*sini)**3/a**3)**Decimal(0.5))-m2

Tsun = 4.925490947e-6
secperday = 3600*24

#data = np.load('1950.before.npz')
data = np.load('1713.all.npz')

fitpars = list(data['fitpars'])
res = data['res']

plist = fitpars
MarkovChain = res[:,:-1]
MCMCSize, parSize= MarkovChain.shape
print 'MCMCSize', MCMCSize

#psr = libstempo.tempopulsar(parfile='before_DDGR.par', timfile='J1950+2414_T2.tim') 
#psr = libstempo.tempopulsar(parfile='after_DDGR.par', timfile='1950.all.tim')
#psr.fit()
#pars = psr.pars() 
#vals = psr.vals()
#errs = psr.errs()
#fitidx = [i for i,p in enumerate(pars) if p in fitpars ]
#vals0 = vals[fitidx]
#errs0 = errs[fitidx]
#parsize = vals0.size

#def extract(par):
    #ipar = plist.index(par)
    #return MarkovChain[:,ipar]*errs0[ipar] + vals0[ipar]

def extract(par):
    if par in plist:
        ipar = plist.index(par)
        return MarkovChain[:,ipar]
    else:
        print par
        return np.float(md.__dict__[par][0])

f, ax = plt.subplots()

M2 = extract('M2').astype(float)
KIN = extract('KIN').astype(float)
I = KIN/180.*np.pi
SINI = np.sin(I)
Pb = extract('PB').astype(float)
A1 = extract('A1').astype(float)

M1 = secperday*Pb/2/np.pi*(np.sqrt(Tsun*(M2*SINI)**3/A1**3))-M2

makesubplot2d(ax, M2, M1) 
#ax.set_xlim(0., 0.5)
#ax.set_ylim(0., 3.0)
show()
