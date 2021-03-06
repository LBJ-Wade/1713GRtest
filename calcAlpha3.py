import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from cmath import *
import os,sys

#from psrpbdot import M1
from datatools.tempo import *
from astropy import coordinates as coord
from astropy import constants as const
from tools.Coordinate import *
import numpy.linalg as linalg
from GalacticGeometry import *

import libstempo 
from threadit import spamit

#from Arrow3D import Arrow3D
secperday = 24*3600
secperyear = secperday*365.24218967
#solardist = 8.34 # old
solardist = 8.34 # Reid et al. 2014, (rmb+14)
PI = np.pi
Tsun = np.float128(4.925490947e-6)
R0 = solardist
#AU =  1.469e13
#Msun = 1.9882e33
#G = 6.673e-8
#c = 2.99792458e10
c = const.c.cgs.value
kpc = const.kpc.cgs.value
AU = const.au.cgs.value #1.469e13
Msun = const.M_sun.cgs.value #1.9882e33
G = const.G.cgs.value

from optparse import OptionParser
usage = "usage: %prog [options] arg"
parser = OptionParser()
parser.add_option("-f", '--parfile', dest="parfile", help="par file")
#parser.add_option("-t", '--timfile', dest="timfile", help="tim file")
parser.add_option("-d", '--npzfile', dest="npzfile", help="npz file")
(options, args) = parser.parse_args(args=sys.argv[1:])
print options
parfile = options.parfile
npzfile = options.npzfile
#timfile = options.timfile
"""load the parfile"""
pf = PARfile(parfile)

"""read some information from the parfile"""
gl, gb = getGpos(pf)
try:
    D = 1./float(pf.PX[0])
except:
    D = float(pf.Dist[0])

if pf.__dict__['BINARY'] in ['DD']:
    E = float(pf.E[0])
    Eerr = float(pf.E[1])
    OMerr = float(pf.OM[1])/180.*np.pi
if pf.__dict__['BINARY'] in ['T2']:
    E = float(pf.ECC[0])
    Eerr = float(pf.ECC[1])
    OMerr = float(pf.OM[1])/180.*np.pi
elif pf.__dict__['BINARY'] in ['ELL1', 'ELL1+', 'T2+']:
    E = np.sqrt(float(pf.EPS1[0]**2 + pf.EPS2[0]**2))
    Eerr = np.sqrt(float(pf.EPS1[1]**2 + pf.EPS2[1]**2))
    e1 = float(pf.EPS1[0])
    e2 = float(pf.EPS2[0])
    er1 = pf.EPS1[1]/pf.EPS1[0]
    er2 = pf.EPS2[1]/pf.EPS2[0]
    er = np.sqrt(float(er1**2 + er2**2))
    OM = np.arctan2(e1,e2)*180./np.pi % 360
    x = e1/e2
    OMerr = 1./(1. + x**2) * er * x * 180./np.pi

"""Galactic acceleration for low z """
#Kz = lambda z:(2.27*z + 3.68*(1-np.exp(-4.31*np.abs(z))) ) * 1.e-9 #Galactic acceleration in z direction (cm/s^2)
def Kz(z):
    absz = np.abs(z)
    sign = z/absz
    return sign * (2.27*absz + 3.68*(1-np.exp(-4.31*absz)) ) * 1.e-9

"""Calculate the coordinate transfermation matrix"""
T, GT = GetTransferMatrix(pf)#, paascnode)

def EdotF(f0, x, px, sini, paascnode, Pb, PMRA, PMDEC, omdot_GR, e1, e2, e1dot, e2dot, m1, m2, w, alpha_1):
    """calculate delta using e1dot and e2dot
    """
    Mtot = m1 + m2
    nb = 2. * np.pi / Pb
    VO = (G * (m1+m2) * Msun * nb) ** (1./3.)

    wserr = 0.9e5 #Kogut et al. 1993, Fixsen et al 1996, Hinshaw et al. 2009
    wsolar = 369.e5 + np.random.randn()*wserr#See ref below (aaa+13 Planck Team: Aghanim, N. et al. 2013. Planck confirms this using a different method)
    lws, bws = 263.99/180.*np.pi, 48.26/180.*np.pi
    w_s = wsolar * (GT.I * np.matrix((np.cos(bws)*np.cos(lws),np.cos(bws)*np.sin(lws),np.sin(bws))).T)
    ws_NSEW = T * w_s
    D = kpc/px
    
    wx = w * (np.matrix((-1., 0., 0.)).T)
    wy = PMRA*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,-1.,0.)).T)
    wz = PMDEC*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,0.,1.)).T)

    w_ns = (wx + wy + wz)
    w =  w_ns + ws_NSEW

    incang = np.arcsin(sini)
    Omgang = paascnode/180.*np.pi
    A = 1. / np.sin(Omgang) / np.tan(incang)
    B = -1./np.tan(Omgang)
    C = -1.

    #C_psr = 0.21*m1 #compactness
    #C_psr = 0.11 * 2 #compactness = sensitivity x 2
    s_p = 0.1

    n_orb = np.matrix((A, B, C)).T
    n_orb= n_orb/linalg.norm(n_orb)
    #print 'n_orb', n_orb.T
    #w_proj = w - n_orb * (w.T*n_orb) 
    w_orb = w - n_orb * (n_orb.T * w)
    #w_proj = w - n_orb * np.dot(w,n_orb) 
    w_leg = linalg.norm(w_orb)
    #w_dir = w_proj/w_leg
    w_dir = w_orb/w_leg
    n_orb = np.array((A, B, C))
    n_orb= n_orb/linalg.norm(n_orb)

    #A_ref = np.matrix((0, -1.* np.sin(Omgang), np.cos(Omgang))) 
    A_ref = np.array((0, -1.* np.sin(Omgang), np.cos(Omgang))) #direction of the accending node, e2, and e2dot
    B_ref = np.cross(n_orb, A_ref) #direction of e1, e1dot

    edot_obs = e1dot * B_ref + e2dot * A_ref
    #e_dir = np.cos(om/180.*np.pi)*A_ref + np.sin(om/180.*np.pi)*B_ref #direction of the eccentrcity e
    e_arr = e1 * B_ref + e2 * A_ref #e1 is the e*sin(omega) part and e2 is the e*cos(omege) part
    ecc = linalg.norm(e_arr)
    edot_GR = omdot_GR * np.cross(n_orb, e_arr)

    edot_a1 = alpha_1 / 4. / c**2  * (m1-m2)/(m1+m2) * nb * VO * w_orb

    edot = edot_obs - edot_GR - edot_a1

    #edot_a3 = alpha_3 * np.pi * s_p * f0 / VO * w_orb

    w_frc = -1. * np.pi * s_p * f0 / VO * w_orb
    #w_frc = 0.25 * C_psr * f0 * Pb * np.cross(np.cross(w_proj, n_orb), n_orb) * np.sqrt(1. - ecc**2) * sini / x / c
    w_frc_norm = linalg.norm(w_frc)
    w_frc_dir = w_frc/w_frc_norm

    return  (np.dot(edot, w_frc_dir)/w_frc_norm)[0][0] #return edot/w_{_|_}


def calcalpha3(F0, A1, PX,  SINI, PAASCNODE, PMRA, PMDEC, M1, M2, PB, ECC, OM, E1DOT, E2DOT, Wr, alpha_1):
    VI = np.arange(M1.size)[np.logical_and(M1 > 1.0, M1 < 2.5)] #valid indices

    OMDOT_GR = 3.*(2*np.pi/PB)**(5./3)*(Tsun*(M1+M2))**(2./3)/(1. - ECC**2)

    alpha3 = []
    for i,sini in enumerate(SINI):
        #print 'i, z[i]', i, z[i], gb, D[i]
        edforce = EdotF(F0, A1[i], PX[i], sini, PAASCNODE[i], PB[i], PMRA, PMDEC, OMDOT_GR[i], E1[i], E2[i], E1DOT[i], E2DOT[i], M1[i], M2[i], Wr, alpha_1[i])
        #edforce = EdotF(F0[i], A1[i], PX[i], sini, PAASCNODE[i], PB[i], OMDOT_GR[i], E1[i], E2[i], E1DOT[i], E2DOT[i], M1[i], M2[i], Wr)
        alpha3.append(edforce)

    return np.array(alpha3)


data = np.load(npzfile)

fitpars = list(data['fitpars'])
res = data['res']

#psr = libstempo.tempopulsar(parfile=parfile, timfile=timfile) 
#psr.fit()
#pars = psr.pars() 
#vals = psr.vals()
#errs = psr.errs()
#fitidx = [i for i,p in enumerate(pars) if p in fitpars ]
#vals0 = vals[fitidx]
#errs0 = errs[fitidx]
#parsize = vals0.size

md = model(parfile)

plist = fitpars
MarkovChain = res[:,:-1]
MCMCSize, parSize= MarkovChain.shape

#def getpar(par):
    #ipar = plist.index(par)
    #return MarkovChain[:,ipar]*errs0[ipar] + vals0[ipar]

def getpar(par):
    if par in plist:
        ipar = plist.index(par)
        return MarkovChain[:,ipar]
    else:
        return np.float(md.__dict__[par][0])

TS99 = lambda pb, p: (pb/p[1])**(1./p[0]) + p[2]
def PbToM2(pb): #based on different theoretical models in Tauris & Savonije 1999
    pars = [(4.5, 1.2e5, 0.12), (4.75, 1.1e5, 0.115), (5.0, 1.e5, 0.11)]
    M2s = np.array([TS99(pb, p) for p in pars])
    return np.random.uniform(M2s.min(), M2s.max())
    
PB = getpar('PB')*secperday  #np.array([float(p[ipb])*secperday for p in MarkovChain])

if "F0" in plist:
    F0 = getpar('F0')
else:
    F0 = np.float(md.F0[0])


if 'M2' in plist:
    M2 = getpar('M2')
else:
    M2 = np.array([PbToM2(p/secperday) for p in PB])

if 'SINI' in plist:
    SINI = getpar('SINI')
elif 'KIN' in plist:
    SINI = np.sin(getpar('KIN')/180.*np.pi)
else:
    COSI = np.random.uniform(0., 1., MCMCSize)
    SINI = np.sqrt(1. - COSI**2)

a = getpar('A1')

if 'E' in plist:
    ECC = getpar('E')
elif 'ECC' in plist:
    ECC = getpar('ECC')
elif 'EPS1' in plist:
    E1 = getpar('EPS1')
    E2 = getpar('EPS2')
    ECC = np.sqrt(E1**2 + E2**2)

try:
    PX = getpar('PX')
except:
    PX = 1./np.array(float(pf.Dist[0]) + np.random.rand(MCMCSize)*float(pf.Dist[1]))

if 'PAASCNODE' in plist:
    PAASCNODE = getpar('PAASCNODE')
elif 'KOM' in plist:
    PAASCNODE = getpar('KOM')
else:
    PAASCNODE = np.random.uniform(0., 360.,  MCMCSize)

if pf.__dict__['BINARY'] in ['DD', 'T2']:
    OM = getpar('OM')
elif pf.__dict__['BINARY'] in ['ELL1', 'ELL1+', 'T2+']:
    E1 = getpar('EPS1')
    E2 = getpar('EPS2')
    OM = np.arctan2(E1,E2)*180./np.pi % 360

if 'EPS1DOT' in plist:
    E1DOT = getpar('EPS1DOT')
if 'EPS2DOT' in plist:
    E2DOT = getpar('EPS2DOT')

PMRA = getpar('PMRA')
PMDEC = getpar('PMDEC')

#M1 = (PB/2/pi*np.sqrt(G*(M2*SINI)**3/a**3)-M2)/Msun
M1 = PB/2/pi*(np.sqrt(Tsun*(M2*SINI)**3/a**3))-M2
#print 'M1', M1

#print 'EPS1DOT:', E1DOT.mean() , E1DOT.std()
#print 'parfile:', vals0[plist.index('EPS1DOT')], errs0[plist.index('EPS1DOT')]
#print 'EPS2DOT:', E2DOT.mean() , E2DOT.std()
#print 'parfile:', vals0[plist.index('EPS2DOT')], errs0[plist.index('EPS2DOT')]

#OMDOT_GR = 3.*(2*np.pi/PB)**(5./3)*(Tsun*(M1+M2))**(2./3)/(1. - ECC**2)
#incang = np.arcsin(SINI.mean())
#Omgang = PAASCNODE.mean()/180.*np.pi
#A = -1./np.tan(incang)
#B = -1./np.tan(Omgang)
#C = 1.
#n_orb = np.array((A, B, C))
#n_orb= n_orb/linalg.norm(n_orb)
#A_ref = np.array((0, -1.* np.sin(Omgang), np.cos(Omgang))) #direction of the accending node, e2, and e2dot
#B_ref = np.cross(n_orb, A_ref) #direction of e1, e1dot
#e1, e2 = E1.mean(), E2.mean()
#e1dot, e2dot = E1DOT.mean(), E2DOT.mean()
#edot_obs = e1dot*B_ref + e2dot*A_ref
#e_arr = e1 * B_ref + e2 * A_ref #e1 is the e*sin(omega) part and e2 is the e*cos(omege) part
#edot_GR = OMDOT_GR.mean() * np.cross(n_orb, e_arr)
#edot = edot_obs - edot_GR
#print 'edot_obs:', edot_obs, 'edot_GR:', edot_GR, 'edot_exc:', edot

#alpha_1 = 1.7e-5 # 0.4^{+3.7}_{-3.1}e-5
#alpha_1 = np.random.randn(SINI.size) * 1.7e-5 + 0.4e-5
alpha_1 = np.zeros(SINI.size)

Wr = np.linspace(-3000.e5, 3000.e5, num=100)
X = []
Y = []
def calc(w):
    return  calcalpha3(F0, a, PX, SINI, PAASCNODE, PMRA, PMDEC, M1, M2, PB, ECC, OM, E1DOT, E2DOT, w, alpha_1).flatten()


#print calc(0.).shape
#sys.exit(0)

X = spamit(calc, [[w] for w in Wr])

for w in Wr:
    Y.append(w*np.ones(X[0].size))
x = np.hstack(X).astype(float)
y = np.hstack(Y).astype(float)
hist2d(x, y, bins=10)
#scatter(x,y, color='k')
colorbar()
xlabel(r"$\^{\alpha}_3$")
ylabel(r'$v_r$ (cm/s)')
show()

print x.shape, y.shape
#np.savez('alpha3_alpha1', Wr=Wr, Alpha3=X)
np.savez('alpha3_only', Wr=Wr, Alpha3=X)

#alpha3 = calcalpha3(F0, a, PX, SINI, PAASCNODE, PMRA, PMDEC, M1, M2, PB, ECC, OM, E1DOT, E2DOT, 0)
#hist(alpha3)
#show()

