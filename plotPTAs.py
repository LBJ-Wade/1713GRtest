from datatools.tempo import *
from pylab import *

md = model('1713.cut.par')
tf = TOAfile('1713.cut.tim')

Legacy = []
NG = []
EP = [] 
for i,toa in enumerate(tf.toalist):
    if not 'pta' in toa.flags:
        Legacy.append(i)
    elif toa.flags['pta'] == 'NANOGrav':
        NG.append(i)
    elif toa.flags['pta'] == 'EPTA':
        EP.append(i)

fig, (ax1, ax2) = subplots(2, 1, sharex=True)

md.tempo2fit(tf)
md.average(groups = {'Legacy':Legacy, 'NANOGrav':NG, 'EPTA':EP}, lapse = 14)

md.plot('mjd', 'DMX', ax=ax1)
md.plot('mjd', 'averes', LegendOn=True, ax=ax2)
show()

#md.plot('freq', 'res', LegendOn=True)
#show()
