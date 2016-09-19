from datatools.tempo import *
from pylab import *
import matplotlib.pyplot  as plt

#m = model("1713.May.noM4O.par")
#m = model("1713.Sep.T1.par")
#m = model("1713_21yr_simple.par")
#m = model("1713.Jun.OMDOT.par")
#t = TOAfile('1713.May.noM4O.tim')
#t = TOAfile('1713.Jun.tim')
#m = model("Feb.T1.RN.par")
#t = TOAfile('1713.Feb.T2.tim')
m = model("1713.all.par")
t = TOAfile('1713.all.tim')
m.tempo2fit(t)
#m.average()
#m.plot('date', 'averes', LegendOn=True)
#show()
oldgroups = m.groups
newgroups = {}
newgroups['ASP-L'] = m.groups['ASP-L']# + m.groups['ASP-S']
newgroups['GASP-8'] = m.groups['GASP-8'] #+ m.groups['GASP-L']
newgroups['Mark3'] = m.groups['M3A-L'] + m.groups['M3B-L']
newgroups['Mark4-L'] = m.groups['M4-L'] #+ m.groups['M42380']
newgroups['Mark4O-L'] = m.groups['M4O-L']
newgroups['ABPP-L'] = m.groups['ABPP-L'] #+ m.groups['ABPP2380']
newgroups['GUPPI-8'] = m.groups['GUPPI-8'] #+m.groups['GUPPI-L']
newgroups['PUPPI-L'] = m.groups['PUPPI-L'] #+m.groups['PUPPI-S']
colorgrp = {'Mark3':'c', 'Mark4-L':'m','ABPP-L':'y', 'ASP-L':'g','GASP-8':'b', 'GUPPI-8':'r', 'PUPPI-L':'k', 'Mark4O-L':'r'}
oldgroups = m.groups
m.groups = newgroups
m.average(groups = newgroups)
#ax1 = subplot(211)
f, (ax1,ax2) = plt.subplots(2, sharex=True)
m.plot('date', 'averes', colors=colorgrp, LegendOn=True, ax=ax1)
ax1.legend(loc=2, numpoints=1, bbox_to_anchor=(1.02,0,1.4,1), ncol=1, borderaxespad=0)
#ax1.legend(bbox_to_anchor=(0.0, 1.02,1, 1.02), numpoints=1, ncol=7, mode='expand', borderaxespad=0)
print m.avewrms
m.groups=oldgroups
newgroups = {}
newgroups['ASP-S'] = m.groups['ASP-S']# + m.groups['ASP-L']
newgroups['GASP-L'] = m.groups['GASP-L'] #+ m.groups['GASP-800']
newgroups['Mark3'] = m.groups['M3A-L'] + m.groups['M3B-L']
newgroups['Mark4-S'] = m.groups['M4-S'] #+ m.groups['M41410']
newgroups['Mark4O-S'] = m.groups['M4O-S']
newgroups['ABPP-S'] = m.groups['ABPP-S'] #+ m.groups['ABPP1410']
newgroups['GUPPI-L'] = m.groups['GUPPI-L']#+m.groups['GUPPI-800']
newgroups['PUPPI-S'] = m.groups['PUPPI-S']#+m.groups['PUPPI-L']
colorgrp = {'Mark3':'c', 'Mark4-S':'m','ABPP-S':'y', 'ASP-S':'g','GASP-L':'b', 'GUPPI-L':'r', 'PUPPI-S':'k', 'Mark4O-S':'r'}
#ax2 = subplot(212)
m.groups = newgroups
m.average()
m.plot('date', 'averes', colors=colorgrp, ax=ax2)
ax2.legend(loc=2, numpoints=1, bbox_to_anchor=(1.02,0, 1.4,1), ncol=1, borderaxespad=0)
#ax2.legend(bbox_to_anchor=(0.0, -.02,1, -.02), numpoints=1, ncol=7, mode='expand', borderaxespad=0)
ax1.set_ylim((-3,3))
ax2.set_ylim((-3,3))
f.subplots_adjust(hspace=0.0)
show()
print m.chisq, m.dof
