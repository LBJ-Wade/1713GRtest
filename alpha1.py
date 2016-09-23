import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
from datatools.tempo import *
from astropy import coordinates as coord
from astropy import constants as const
from tools.Coordinate import *
from psrpbdot import M1
import numpy.random as npr
from cmath import *
import os,sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy.linalg as linalg
from GalacticGeometry import *

PI = np.pi
Tsun = 4.925490947e-6
secperday = 24*3600
secperyear = secperday*365.24218967
c = const.c.cgs.value
kpc = const.kpc.cgs.value
AU = const.au.cgs.value #1.469e13
Msun = const.M_sun.cgs.value #1.9882e33
G = const.G.cgs.value


"""
Shao and Wex 2012
"""
alpha_1 = 1.7e-5 # 0.4^{+3.7}_{-3.1}e-5
alpha_2 = 1.8e-4 # 95% confidence limit
alpha_3 = 1.e-21 # benchmark
s_p = 0.1

psr = model('1713.cut.par')

T, GT = GetTransferMatrix(psr)

m1 = float(M1(psr))
m2 = float(psr.M2[0])
pb = float(psr.PB[0]) * secperday
ecc = np.sqrt(float(psr.EPS1[0])**2 + float(psr.EPS2[0])**2)
spin = float(psr.F0[0]) 

q = m1/m2

print 'ecc:', ecc, 'q:', q

nb = 2. * PI / pb

omdot_GR = 3.*nb**(5./3)*(Tsun*(m1+m2))**(2./3)/(1. - ecc**2)
print 'omdot_GR:', omdot_GR * secperyear * 57.2958

Fe = 1./(1. + np.sqrt(1-ecc**2))

VO = (G * (m1+m2) * Msun * nb) ** (1./3.)
print 'VO:', VO

#w = 1000.e5
w = -1000.e5


wserr = 0.9e5 #Kogut et al. 1993, Fixsen et al 1996, Hinshaw et al. 2009
wsolar = 369.e5 + np.random.randn()*wserr#See ref below (aaa+13 Planck Team: Aghanim, N. et al. 2013. Planck confirms this using a different method)
lws, bws = 263.99/180.*np.pi, 48.26/180.*np.pi
w_s = wsolar * (GT.I * np.matrix((np.cos(bws)*np.cos(lws),np.cos(bws)*np.sin(lws),np.sin(bws))).T)
ws_NSEW = T * w_s

px = float(psr.PX[0])
D = kpc/px
PMRA = float(psr.PMRA[0])
PMDEC = float(psr.PMDEC[0])
wx = w * (np.matrix((-1., 0., 0.)).T)
wy = PMRA*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,-1.,0.)).T)
wz = PMDEC*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,0.,1.)).T)
w_ns = (wx + wy + wz)
w =  w_ns + ws_NSEW

incang = float(psr.KIN[0])/180.*np.pi
Omgang = float(psr.KOM[0])/180.*np.pi
A = 1. / np.sin(Omgang) / np.tan(incang)
B = -1./np.tan(Omgang)
C = -1.

n_orb = np.matrix((A, B, C)).T
n_orb = n_orb/linalg.norm(n_orb)
print 'n_orb:', [x for x in n_orb]
print 'w_solar:', [x for x in ws_NSEW]
print 'w_psr:', [x for x in w_ns]
print 'w = w_psr + w_solar:', w.flatten()
print 'w - (w dot n_orb)*n_orb:', (w - n_orb * (n_orb.T * w)).flatten()

w_orb = w - n_orb * (n_orb.T * w)
edot_a1 = alpha_1 / 4. / c**2  * (m1-m2)/(m1+m2) * nb * VO * w_orb
edot_a3 = alpha_3 * PI * s_p * spin / VO * w_orb
#print edot_a1, edot_a3

e_x = np.matrix((0., -1.*np.sin(Omgang), np.cos(Omgang))).T
e_y = np.matrix(np.cross(np.array(n_orb).flatten(), np.array(e_x).flatten())).T

print 'e_x_dir, e_y_dir', e_x.flatten(), e_y.flatten()
print 'edot_a1:', 'EPSDOT1', edot_a1.T * e_x , 'EPSDOT2', edot_a1.T * e_y
print 'edot_a3:', 'EPSDOT1',edot_a3.T * e_x , 'EPSDOT2', edot_a3.T * e_y

