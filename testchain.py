import numpy as np

data = np.genfromtxt('chain_1.txt', dtype='f')
pars = np.genfromtxt('pars.txt', dtype='|S30')
#print pars.size, pars.shape
#print data.size, data.shape
np.savez('1713.all', fitpars=pars, res=data[...,:59])

#npzfile = '1713.all.npz'
#data = np.load(npzfile)

#fitpars = list(data['fitpars'])
#res = data['res']

#print len(fitpars), res.shape

#m2idx = fitpars.index('PB')
#print res[:,m2idx]
