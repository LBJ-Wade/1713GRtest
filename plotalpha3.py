from pylab import *
import matplotlib.pyplot as plt
from datatools.tempo import initmatplotlib
from datatools.tempo import *
from astropy import constants as const

import numpy as np
import scipy.stats as ss

mpl = initmatplotlib(cols=0)

kpc = const.kpc.cgs.value
secperday = 24*3600
secperyear = secperday*365.24218967

pf = PARfile('1713.cut.par')
PMRA = float(str(pf.PMRA[0]))
PMDEC = float(str(pf.PMDEC[0]))
px = float(str(pf.PX[0]))
D = kpc/px

wy = PMRA*1.e-3/60./60.*np.pi/180./secperyear * D
wz = PMDEC*1.e-3/60./60.*np.pi/180./secperyear * D 
wpm = np.sqrt(wy**2 + wz**2)

data = np.load('alpha3_result.npz')
Wr = data['Wr']
Alpha3 = data['Alpha3']

size = Alpha3[0].size
idx = np.arange(size)
#lowlim1, hilim1 = int(0.16*size), int(0.84*size)
#lowlim2, hilim2 = int(0.025*size), int(0.975*size)
lowlim1, hilim1 = int(0.025*size), int(0.975*size)
lowlim2, hilim2 = int(0.0015*size), int(0.9985*size)
braket1 = []
braket2 = []
normbk1 = []
normbk2 = []
for i,w in enumerate(Wr):
    a3 = Alpha3[i]
    a3.sort()
    #print a3[lowlim1], a3[hilim1]
    braket1.append((a3[lowlim1], a3[hilim1]))
    #print a3[lowlim1], a3[hilim1]
    braket2.append((a3[lowlim2], a3[hilim2]))
    mu, std = ss.norm.fit(a3)
    #print mu, std
    #print ss.norm.interval(0.95, loc=mu, scale=std)
    #normbk1.append(ss.norm.interval(0.68, loc=mu, scale=std))
    normbk1.append(ss.norm.interval(0.95, loc=mu, scale=std))
    normbk2.append(ss.norm.interval(0.997, loc=mu, scale=std))

braket1 = np.array(braket1)
braket2 = np.array(braket2)
normbk1 = np.array(braket1)
normbk2 = np.array(braket2)

#plot(braket1[:,0], Wr, 'r-')
#plot(braket1[:,1], Wr, 'r-')
#plot(braket2[:,0], Wr, 'r--')
#plot(braket2[:,1], Wr, 'r--')

print '95% limits:', normbk1[:,0].min(), normbk1[:,1].max()
#print wpm/1.e8

plot(normbk1[:,0]/1.e-21, Wr/1.e8, 'k-', linewidth=2)
plot(normbk1[:,1]/1.e-21, Wr/1.e8, 'k-', linewidth=2)
plot(normbk2[:,0]/1.e-21, Wr/1.e8, 'k--', linewidth=2)
plot(normbk2[:,1]/1.e-21, Wr/1.e8, 'k--', linewidth=2)
#plt.axhline(y=wpm/1.e8, linewidth=2)
#plt.axhline(y=wpm/-1.e8, linewidth=2)
xlabel(r"$\^{\alpha}_3 $ $(10^{-21})$")
ylabel(r'$v_r$ ($10^{8}$cm/s)')
show()
