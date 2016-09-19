import matplotlib.pylab as plt
import numpy as np
from datatools.tempo import *
from numpy import genfromtxt
from datatools.MJD import MJD_to_datetime
from pylab import *



#m = model('1713.Apr.dmx.par')
#t = TOAfile('1713.Apr.tim')
#m = model('Feb.T1.RN.par')
#t = TOAfile('1713.Feb.T2.tim')
m = model('1713.all.par')
t = TOAfile('1713.all.tim')
m.tempo2fit(t)
#m.tempofit(t) #, GLS=True)
#os.system('tempo -f %s %s -a' % (m.parfile, t.toafile)) #run tempo to generate the chisun.tmp file
#sys.exit(0)
#phisun = genfromtxt('phisun.tmp')
#phisun = phisun[:len(t.toalist)]
DMX, DMXErr, DMXR1, DMXR2 = m.dmxlist
meanDMX = mean([float(DMX[i]) for i in DMX])
dmxvalues = genfromtxt('PaulDMX.par', dtype = [('DMXEP', 'f8'), ('DMX', 'f8'), ('DMXErr', 'f8'), ('DMXR1', 'f8'), ('DMXR2', 'f8'), ('DMXF1', 'f8'), ('DMXF2', 'f8'), ('DMXbin', 'S5')])
DMX= {}
DMXErr = {}
DMXR1 = {}
DMXR2 = {}
DM = float(m.DM)
#print meanDMX - DM
for rec in dmxvalues:
    i = int(rec['DMXbin'][2:].lstrip('0'))
    DMX[i] = rec['DMX'] + meanDMX
    DMXErr[i] = rec['DMXErr']
    DMXR1[i] = rec['DMXR1']
    DMXR2[i] = rec['DMXR2']

def MJD_to_year(mjd):
    d = MJD_to_datetime(mjd)
    return d.year + (float(d.strftime("%j")) + float(d.strftime("%H"))/24)/365.24218967

DMXgrp = {}
for j in DMXR1.keys():
    DMXgrp[j] = []
for i in range(len(t.toalist)):
    for j in DMXR1.keys():
        if float(t.toalist[i].TOA) > float(DMXR1[j]) and float(t.toalist[i].TOA) < float(DMXR2[j]):
            DMXgrp[j].append(i)

DMXR = {}
for j in DMXR1.keys():
    if not j == 1:
        DMXR[j] = (DMXR1[j] + DMXR2[j])/2
DMidx = []
for j in DMXR:
    if not j == 1:
        DMidx.append((DMXR[j], j))
DMidx = np.array(DMidx, dtype=[('dm', float ), ('i', int)] )
DMidx.sort(order='dm')
idx = DMidx['i']

for j in DMXR1.keys():
    #DMXR[j] = (DMXR1[j] + DMXR2[j])/2
    toa_in_DMbin = []
    toa_sigma = []
    for i in range(len(m.toa)):
        if m.toa[i] > float(DMXR1[j]) and m.toa[i] < float(DMXR2[j]):
            toa_in_DMbin.append(m.toa[i])
            toa_sigma.append(m.err[i])
    if len(toa_in_DMbin) > 0:
        toa_in_DMbin = np.array(toa_in_DMbin)
        toa_sigma = np.array(toa_sigma)
        weights = 1/toa_sigma**2
        weighted_average = sum(weights*toa_in_DMbin)/sum(weights)
        DMXR[j] = weighted_average 
        #print len(toa_in_DMbin), weighted_average 
    else:
        DMXR[j] = (DMXR1[j] + DMXR2[j])/2



#dmxphi = np.array([np.mean([phisun[i] for i in DMXgrp[j]]) for j in DMXgrp.keys() if DMXR1[j]>53200 and not len(DMXgrp[j]) == 0 and (DMXR2[j]<54760 or DMXR1[j]>54870)])
#phi = np.array([np.mean([phisun[i] for i in DMXgrp[j]]) for j in DMXgrp.keys() if DMXR1[j]>53200 and not len(DMXgrp[j]) == 0])
#phi = np.array([np.mean([phisun[i] for i in DMXgrp[j]]) for j in DMXgrp.keys() if not len(DMXgrp[j]) == 0])
#phiidx = [ j for j in DMXgrp.keys() if DMXR1[j]>53200 and not len(DMXgrp[j]) == 0 ]

#phiidx = [ j for j in DMXgrp.keys() if not len(DMXgrp[j]) == 0 ]

#phidmr = np.array([float(DMXR[i])-50000 for i in phiidx]) 
#data = np.vstack((phidmr, phi)).T.flatten()
#dataview = data.view(dtype=[('dmr', np.float), ('phi', np.float)])
#dataview.sort(order='dmr')
#phidmr = dataview['dmr']
#phi = dataview['phi']

dmr = np.array([float(DMXR[i])-50000 for i in idx]) 
dmx = np.array([float(DMX[i])-16. for i in idx]) 
dmxerr = np.array([float(DMXErr[i]) for i in idx])
dmr1 = dmr[dmr<3200]
dmr2 = dmr[dmr>3200]
dmx1 = dmx[dmr<3200]
dmx2 = dmx[dmr>3200]
dmxerr1 = dmxerr[dmr<3200]
dmxerr2 = dmxerr[dmr>3200]

#f, (ax1, ax2) = plt.subplots(2,1, sharex=True)
f, ax = plt.subplots()

ax.errorbar([MJD_to_year(d+50000) for d in dmr1], dmx1, dmxerr1, fmt='k.')
#ax1.set_ylim(15.9688, 15.9705)

ax.errorbar([MJD_to_year(d+50000) for d in dmr2], dmx2, dmxerr2, fmt='k.')
#ax2.set_ylim(15.9626, 15.9639)
#ax.set_ylim(-0.0338, -0.0315)
ax.set_ylim(-0.0308, -0.0285)

#ax1.spines['bottom'].set_visible(False)
#ax2.spines['top'].set_visible(False)
#ax1.xaxis.tick_top()
#ax1.tick_params(labeltop='off')
#ax2.xaxis.tick_bottom()
#ax1.set_ylim((-0.0034, -0.0019))
#ax2.set_ylim((-0.0096, -0.0081))
#ax1.set_xlim((0, 7000))
#ax2.set_xlim((0, 7000))

#d= .015
#kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
#ax1.plot((-d,+d), (-d,+d), **kwargs)
#ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
#kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False )
#ax2.plot((-d,+d), (1-d,1+d), **kwargs)
#ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

#ax.annotate('Instrumental', fontsize=18, xy=(MJD_to_year(53200), -0.0087), xycoords='data', xytext=(-59,-90), textcoords='offset points', arrowprops=dict(arrowstyle="->"))

#plt.ylabel('                                                DM - 16 (pc~mc$^{-3}$)')
plt.ylabel('DM - 16 (pc~cm$^{-3}$)')
#plt.xlabel('MJD-50000')
plt.xlabel('Year')
#show()
#sys.exit()

"""plot solar elongation """
import ephem
from MJD import *
from astropy import coordinates as coord
from tools.Coordinate import RA, Dec
mjd = np.linspace(50900, 56500, 500)
#mjd = np.linspace(53200, 57000, 500)
#dates = [MJD_to_datetime(m).strfmt('%Y/%m/%d') for  m in mjd]
dates = [MJD_to_datetime(jd) for  jd in mjd]
if type(m.RAJ) == list:
    ra = RA(m.RAJ[0])
    dec = Dec(m.DECJ[0])
    #pulsar = ephem.Equatorial(m.RAJ[0], m.DECJ[0], epoch='2000')
    #psr = coord.FK5Coordinates(str(ra) +' '+str(dec))
    #print str(ra) +' '+str(dec)
    psr = coord.SkyCoord(str(ra) +' '+str(dec))
PHI = []
pr = psr.ra.rad
pd = psr.dec.rad
for date in dates:
    sun = ephem.Sun()
    sun.compute(date, epoch='2000')
    sunstr = str(RA(str(sun.g_ra)))+' '+ str(Dec(str(sun.g_dec)))
    sunpos = coord.SkyCoord(sunstr)
    sr = sunpos.ra.rad
    sd = sunpos.dec.rad
    psr = np.array((np.cos(pd)*np.sin(pr), np.cos(pd)*np.cos(pr), np.sin(pd)))
    sun = np.array((np.cos(sd)*np.sin(sr), np.cos(sd)*np.cos(sr), np.sin(sd)))
    PHI.append(np.arccos(np.dot(psr, sun))*180./np.pi)
PHI = np.array(PHI)
#ax2rax = ax2.twinx()
ax2rax = ax.twinx()
ax2rax.spines["right"].set_visible(True)
#ax2rax.plot(phidmr, phi, 'y:', label="Solar Elongation")
ax2rax.plot([MJD_to_year(d) for d in mjd], PHI, 'k:', label="Solar Elongation")
#ax2rax.set_xlim((0, 7000))
ax2rax.set_ylabel("Solar Elongation (deg)")
ax2rax.set_ylim(0, 240)
ax2rax.set_yticks([0,  60, 120,  180])

from scipy.optimize import curve_fit , leastsq
data = np.load('strucfunc.npy')
delays = data[0,...]
ddmsqs = data[1,...]
biners = data[2,...]
ddmers = data[3,...]
#A,B,r,p,s = linregress(np.log(delays), np.log(ddmsqs))
a = axes([.62, .6, .25, .3])
def func(x,  b , c):
    return c * (x)**(b-2) 
popt, pcov = curve_fit(func, delays, ddmsqs)
popt, pcov = curve_fit(func, delays, ddmsqs, sigma=ddmers, p0= popt)
B,C = popt
Berr = np.sqrt(pcov[0,0])
Cerr = np.sqrt(pcov[1,1])
from round import shortform as SF
print 'beta:', SF((B,Berr)), 'C:', SF((C,Cerr))
x = np.linspace(10, 1.e4, 30)
y = func(x, B, C)
a.plot(x,y, 'k-')

#def Kfunc(x,  c):
    #return c*(x)**(11./3-2) 
#popt, pcov = curve_fit(Kfunc, delays, ddmsqs, sigma=ddmers)
#C = popt[0]
#y = Kfunc(x, C)#Kolmogorov
#y = func(x, 11./3, C)#Kolmogorov
#a.plot(x,y, 'k--')

a.errorbar(delays, ddmsqs, xerr=biners, yerr=ddmers, fmt='ok')
a.semilogx()
a.semilogy()
a.set_ylabel(r'D$_{\phi}$ (rad$^2$)')
a.set_xlabel(r'$\tau$ (days)')
a.set_title('Structure Function')

plt.show()

