from datatools.tempo import *

import matplotlib.pyplot as plt
from pylab import *

initmatplotlib(cols=2)

fig, (ax1, ax2, ax3) = subplots(3, 1, sharex=True)

md = model('1713.pbdot.par')
tf = TOAfile('1713.all.tim')

md.tempo2fit(tf)
md.average(groups='all', lapse=14)

#md.plot('mjd', 'DMX', ax=ax1, color='k')
#md.plot('mjd', 'averes', ax=ax2, NoZeroLine=True)
#md.plot('ophase', 'averes', NoZeroLine=True, color='k')
md.plot('mjd', 'averes', NoZeroLine=True, color='k', ax=ax1)
ax1.set_ylim(-2., 2.)

mjd = md.avetoa['all']
averes1 = md.averes['all']
aveerr1 = md.aveerr['all']

#md.PBDOT = 0
#md.write('1713.nopbdot.par')
del md
md = model('1713.nopbdot.par')
md.tempo2fit(tf)
md.average(groups='all', lapse=14)

md.plot('mjd', 'averes', NoZeroLine=True, color='k', ax=ax2)
ax2.set_ylim(-2., 2.)

averes2 = md.averes['all']
aveerr2 = md.aveerr['all']

ax3.errorbar(mjd, averes2-averes1, yerr=np.sqrt(aveerr1**2 + aveerr2**2), color='k')
ax3.set_ylim(-2.,2.)

show()
