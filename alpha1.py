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
alpha_1 = -3e-4 # 0.4^{+3.7}_{-3.1}e-5
alpha_2 = 1.8e-4 # 95% confidence limit


psr = model('1713.cut.par')

m1 = float(M1(psr))
m2 = float(psr.M2[0])
pb = float(psr.PB[0]) * secperday
ecc = np.sqrt(float(psr.EPS1[0])**2 + float(psr.EPS2[0])**2)

q = m1/m2

print 'ecc:', ecc, 'q:', q

nb = 2. * PI / pb

omdot_GR = 3.*nb**(5./3)*(Tsun*(m1+m2))**(2./3)/(1. - ecc**2)
print 'omdot_GR:', omdot_GR * secperyear * 57.2958

Fe = 1./(1. + np.sqrt(1-ecc**2))

VO = (G * (m1+m2) * Msun * nb) ** (1./3.)

print 'VO:', VO

w = 1000.e5
print 'assuming a speed of 1000km/s', w

edot_est = alpha_1 * (q-1)/(q+1) / 4 / c**2 * nb * VO * w
print edot_est

w = 1000.e5
print 'assuming a speed of 1000km/s', w
edot_est = alpha_2 * nb * Fe * w**2 * ecc / c**2
print edot_est
