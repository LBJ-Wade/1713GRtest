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

#from Arrow3D import Arrow3D
secperday = 24*3600
#solardist = 8.34 # old
solardist = 8.34 # Reid et al. 2014, (rmb+14)
PI = np.pi
Tsun = np.float128(4.925490947e-6)
secperday = 24*3600
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

def EdotF(kr, zeta, z, sini, paascnode, omdot_GR, e1, e2, e1dot, e2dot):
    """calculate delta using e1dot and e2dot
    """

    g_r = GT.I * kr * ( np.matrix((np.cos(0.-zeta),np.sin(0.-zeta),0)).T) #X/linalg.norm(X)
    g_z = GT.I * Kz(z) * (np.matrix((0., 0., -1.)).T) 
    g = g_r + g_z
    g_NSEW = T * g

    incang = np.arcsin(sini)
    Omgang = paascnode/180.*np.pi
    A = -1./np.tan(incang)
    B = -1./np.tan(Omgang)
    C = 1.

    n_orb = np.matrix((A, B, C)).T
    n_orb= n_orb/linalg.norm(n_orb)
    #print 'n_orb', n_orb.T

    g_proj = g_NSEW - n_orb * (g_NSEW.T*n_orb) 
    KG  = float(linalg.norm(g_proj))
    g_proj = np.array(g_proj).flatten()
    #A_ref = np.matrix((0, -1.* np.sin(Omgang), np.cos(Omgang))) 
    n_orb = np.array((A, B, C))
    n_orb= n_orb/linalg.norm(n_orb)
    A_ref = np.array((0, -1.* np.sin(Omgang), np.cos(Omgang))) #direction of the accending node, e2, and e2dot
    B_ref = np.cross(n_orb, A_ref) #direction of e1, e1dot
    #print 'n_orb', n_orb, 'A_ref', A_ref, 'g_proj', g_proj
    #g_dir = g_proj/KG
    #g_ang = np.arccos(A_ref * g_dir)
    #g_norm = linalg.norm(g)

    g_frc = np.cross(g_proj, n_orb)
    g_frc_norm = linalg.norm(g_frc)
    g_frc_dir = g_frc/g_frc_norm
    edot_obs = e1dot * B_ref + e2dot * A_ref
    #e_dir = np.cos(om/180.*np.pi)*A_ref + np.sin(om/180.*np.pi)*B_ref #direction of the eccentrcity e
    e_arr = e1 * B_ref + e2 * A_ref #e1 is the e*sin(omega) part and e2 is the e*cos(omege) part
    edot_GR = omdot_GR * np.cross(n_orb, e_arr)
    edot = edot_obs - edot_GR

    return  np.dot(edot, g_frc_dir)/g_frc_norm #return edot/g_{_|_}


def calcdelta(PX, A1, SINI, PAASCNODE, M1, M2, PB, ECC, OM, E1DOT, E2DOT):
    VI = np.arange(M1.size)[np.logical_and(M1 > 1.0, M1 < 2.5)] #valid indices
    #print 'validindices', VI
    
    p = A1 * c * (1 - ECC**2) / SINI
    fac = 1.5 * np.sqrt(p/G/M2/Msun)

    D = 1./PX
    Omega = PAASCNODE/180.*np.pi
    z = np.sin(gb) * D
    R1 = np.sqrt(R0**2 + (D*np.cos(gb))**2 -2 * R0 * D * np.cos(gb) * np.cos(gl))
    coszeta = (R0**2 + R1**2 - D**2 + z**2)/R0/R1/2.
    zeta = np.arccos(coszeta)
    #print 'PX, D, z, gb', PX.mean(), D.mean(), z.mean(), gb
    #print 'Galactic acceleration in z direction in cm/s^2', z, Kz(z)

    Mtot =  M1 + M2
    #Omega_G = 27.2 #km s^-1 kpc^-1
    Omega_G = 30.57 #km s^-1 kpc^-1 +/- 5.1 Ref: Reid et al. 2014 (rmb+14)
    kpcinkm = 3.0857e16
    Kr =  Omega_G**2 * R0**2 / R1 / kpcinkm * 1.e5 #Galactic acceleration in radio direction (cm/s^2)

    OMDOT_GR = 3.*(2*np.pi/PB)**(5./3)*(Tsun*(M1+M2))**(2./3)/(1. - ECC**2)
    #print 'Galactic acceleration in radio direction in cm/s^2', Kr
    #KGarray = []
    #THETA = []
    deltas = []
    for i,sini in enumerate(SINI):
        #print 'i, z[i]', i, z[i], gb, D[i]
        edforce = EdotF(Kr[i], zeta[i], z[i], sini, PAASCNODE[i], OMDOT_GR[i], E1[i], E2[i], E1DOT[i], E2DOT[i])
        deltas.append(edforce/fac[i])
        #THETA.append(theta)
        #KGarray.append(kg)
    #THETA = np.array(THETA)
    #THETA[THETA<0]+=(np.pi*2)
    #KG = np.array(KGarray)
    #EF = Delta * ( 0.5 * KG / Mtot /Tsun/c/(2*PI/PB)**2 )

    return np.array(deltas)


data = np.load(npzfile)

fitpars = list(data['fitpars'])
res = data['res']

#psr = libstempo.tempopulsar(parfile='1713.all.omdot.par', timfile='1713.all.tim') 
#psr.fit()
#pars = psr.pars() 
#vals = psr.vals()
#errs = psr.errs()
#fitidx = [i for i,p in enumerate(pars) if p in fitpars ]
##print fitidx
#vals0 = vals[fitidx]
#errs0 = errs[fitidx]
#parsize = vals0.size

md = model(parfile)

#MarkovChain = res[:,1:]
plist = fitpars
MarkovChain = res[:,:-1]
MCMCSize, parSize= MarkovChain.shape

MCMCSize = len(MarkovChain)

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

deltas = calcdelta(PX, a, SINI, PAASCNODE, M1, M2, PB, ECC, OM, E1DOT, E2DOT)
hist(deltas, 30, normed=1, histtype='step')
xlabel("$\Delta$")
show()


print deltas.shape
np.savez('delta_result', delta=deltas)
