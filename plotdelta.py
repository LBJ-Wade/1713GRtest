from pylab import *
import matplotlib.pyplot as plt
from datatools.tempo import initmatplotlib

import numpy as np
import scipy.stats as ss

mpl = initmatplotlib(cols=1)

data = np.load('delta_result.npz')
delta = data['delta']

size = delta.size
idx = np.arange(size)
lowlim1, hilim1 = int(0.16*size), int(0.84*size)
lowlim2, hilim2 = int(0.025*size), int(0.975*size)
braket1 = []
braket2 = []
normbk1 = []
normbk2 = []
delta.sort()
braket1=(delta[lowlim1], delta[hilim1])
braket2=(delta[lowlim2], delta[hilim2])
mu, std = ss.norm.fit(delta)
#normbk1=ss.norm.interval(0.68, loc=mu, scale=std)
#normbk2=ss.norm.interval(0.95, loc=mu, scale=std)
normbk1=ss.norm.interval(0.95, loc=mu, scale=std)
normbk2=ss.norm.interval(0.997, loc=mu, scale=std)

braket1 = np.array(braket1)
braket2 = np.array(braket2)
normbk1 = np.array(braket1)
normbk2 = np.array(braket2)

frame = plt.gca()
hist(delta, 20, normed=1, histtype='step', color='k', linewidth=2)
vlines(normbk1, 0, 800, color='k', linestyles='solid', linewidth=2)
vlines(normbk2, 0, 800, color='k', linestyles='dashed', linewidth=2)
#plt.axvline(normbk1[0], color='k', linestyles='solid')
#plt.axvline(normbk1[1], color='k', linestyles='solid')
#plt.axvline(normbk2[0], color='k', linestyles='dash')
#plt.axvline(normbk2[1], color='k', linestyles='dash')
xlim(-0.002, 0.002)
ylim(0, 800)
frame.axes.get_yaxis().set_visible(False)
frame.axes.get_xaxis().set_ticks([-0.002, -0.001, 0, 0.001, 0.002])
frame.axes.get_xaxis().set_ticklabels(['-0.002', '-0.001', '0', '0.001', '0.002'])
#tick.set_visible(False)
xlabel("$\Delta$")
show()
