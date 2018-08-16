import numpy as np
from round import TexStyle as SF

#data = np.genfromtxt("g_rob.txt", dtypes=[('gx'), (), (), ()])
data =np.loadtxt("g_orb.txt", dtype='f')

print data.shape

g_x = data[:,0]
g_y = data[:,1]
print (g_x.mean(), np.std(g_x)), (g_y.mean(), np.std(g_y))
print SF([g_x.mean(), np.std(g_x)]), SF([g_y.mean(), np.std(g_y)])
