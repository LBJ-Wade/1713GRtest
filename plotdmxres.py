from datatools.tempo import *

import matplotlib.pyplot as plt
from pylab import *

initmatplotlib(cols=2)

fig, (ax1, ax2) = subplots(2, 1, sharex=True)

md = model('1713.cut.par')
tf = TOAfile('1713.cut.tim')

md.tempo2fit(tf)

md.average()

md.plot('mjd', 'DMX', ax=ax1, color='k')
md.plot('mjd', 'averes', ax=ax2, NoZeroLine=True)

show()
