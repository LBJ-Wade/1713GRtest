from datatools.tempo import *
from TexTable import deluxetable
from datatools.MJD import MJD_to_datetime
from pylab import *
from matplotlib.patches import Rectangle as Box
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib

t = TOAfile('1713.Feb.T2.tim')
m = model('Feb.T1.RN.par')

newgroups = {}
newgroups['Mark III (L)'] = t.groups['M3A-L'] + t.groups['M3B-L']
newgroups['Mark IV (L)'] = t.groups['M4-L']
newgroups['Mark IV (S)'] = t.groups['M4-S']
newgroups['Mark IV-O (L)'] = t.groups['M4O-L']
newgroups['Mark IV-O (S)'] = t.groups['M4O-S']
newgroups['ABPP (L)'] = t.groups['ABPP-L']
newgroups['ABPP (S)'] = t.groups['ABPP-S']
newgroups['ASP (L)'] = t.groups['L-wide_ASP']
newgroups['ASP (S)'] = t.groups['S-wide_ASP']
newgroups['GASP (820)'] = t.groups['Rcvr_800_GASP']
newgroups['GASP (L)'] = t.groups['Rcvr1_2_GASP']
newgroups['GUPPI (820)'] = t.groups['Rcvr_800_GUPPI']
newgroups['GUPPI (L)'] = t.groups['Rcvr1_2_GUPPI']
newgroups['PUPPI (L)'] = t.groups['L-wide_PUPPI']
newgroups['PUPPI (S)'] = t.groups['S-wide_PUPPI']

oldgroups = t.groups
t.groups = newgroups
m.tempofit(t)
m.average(groups=t.groups)
 
keys = ['Mark III (L)', 'Mark IV (L)', 'Mark IV (S)', 'Mark IV-O (L)', 'Mark IV-O (S)', 'ABPP (L)','ABPP (S)','ASP (L)','ASP (S)','GASP (820)', 'GASP (L)', 'GUPPI (820)', 'GUPPI (L)', 'PUPPI (L)', 'PUPPI (S)']

BandWidth = [40, 10, 10, 10, 10, 56, 112, 64, 64, 64, 64, 200, 730, 615, 355]
Integration = [47, 58, 29, 60, 30, 60, 30, 20, 20, 20, 20, 20, 20, 20, 20]

NoOfTOA = []
for k in keys:
    NoOfTOA.append(len(t.groups[k]))

DateRange = []

def datetoyear(dt):
    year = dt.strftime("%Y")
    day = dt.strftime("%j")
    return float(year) + float(day)/365.25

def findrange(grp):
    toas = [t.toalist[i] for i in t.groups[grp]]
    mjds = [float(toa.TOA) for toa in toas]
    maxd = max(mjds)
    mind = min(mjds)
    #return '%s$-$%s' % (MJD_to_datetime(mind).strftime('%Y %b'), MJD_to_datetime(maxd).strftime('%Y %b'))
    #print datetoyear(MJD_to_datetime(mind)), datetoyear(MJD_to_datetime(maxd))
    return datetoyear(MJD_to_datetime(mind)), datetoyear(MJD_to_datetime(maxd))



for k in keys:
    krange = findrange(k)
    #print k, krange
    DateRange.append(krange)
NoOfEpoch = []
for k in keys:
    NoOfEpoch.append(len(m.averes[k]))

CenterFreq = [1400, 1410, 2380, 1410, 2380, 1410, 2380, 1410, 2350, 800, 1410, 820, 1515, 1457, 2227]

#param = matplotlib.rcParams
#color_circle = param['axes.color_cycle']
#color_circle_len = len(color_circle)
#colors = [color_circle[i % color_circle_len] for i in range(len(keys))]

#colors = [float(i+1)/len(keys) for i in range(len(keys))]
colors = []
for i in range(len(keys)):
    if keys[i].startswith("Mark III"):
        colors.append(0.1)
    elif keys[i].startswith("Mark IV"):
        colors.append(0.2)
    elif keys[i].startswith("ABPP"):
        colors.append(0.3)
    elif keys[i].startswith("ASP"):
        colors.append(0.4)
    elif keys[i].startswith("PUPPI"):
        colors.append(0.5)
    elif keys[i].startswith("GASP"):
        colors.append(0.8)
    else: #keys[i].startswith("GUPPI"):
        colors.append(0.9)

        

fig, ax=plt.subplots()
Boxes = []
labeldict = {}
colors = np.array(colors + [colors[-1]])
#CC = matplotlib.cm.ScalableMap()
#CC = plt.get_cmap('jet')
for i,key in enumerate(keys):
    timerange = DateRange[i]
    timemin, timemax = timerange
    bandwid = BandWidth[i]
    cfreq = CenterFreq[i]
    box = Box((timemin, cfreq-bandwid/2), timemax-timemin, bandwid)
    if key.startswith("Mark III"):
        labeldict["Mark III"] = box
    elif key.startswith("Mark IV"): 
        labeldict["Mark IV"] = box
    elif key.startswith("ABPP"): 
        labeldict["ABPP"] = box
    elif key.startswith("ASP"): 
        labeldict["ASP"] = box
    elif key.startswith("GASP"): 
        labeldict["GASP"] = box
    elif key.startswith("PUPPI"): 
        labeldict["PUPPI"] = box
    elif key.startswith("GUPPI"): 
        labeldict["GUPPI"] = box
    Boxes.append(box)

Boxes.append(Box((timemin, 1770), timemax-timemin, 110))
p = PatchCollection(Boxes, cmap=matplotlib.cm.jet, alpha=0.5)
p.set_array(colors)
ax.add_collection(p)

ax.set_xlim((1990, 2015))
ax.set_ylim((300, 2600))
ax.set_xlabel("year")
ax.set_ylabel("frequency (MHz)")

#labels = ("Mark III", "Mark IV", "ABPP", "ASP", "GASP", "PUPPI", "GUPPI")
#legend = plt.legend([labeldict[l] for l in labels], labels,  loc="upper left", labelspacing=0.1, fontsize="small")
#plt.setp(legend.get_texts(), fontsize="small")

plt.show()

#Caption = '21 year J1713+0747 observations.'
#colnames = ['System', 'Dates', 'Number of', 'Epochs', 'Bandwidth', 'Typical Integration']
#data = [keys, DateRange, NoOfTOA, NoOfEpoch, BandWidth, Integration]

#table = deluxetable(Caption=Caption, colnames=colnames, data=data, label='tab:obs', colsetting='lccccc')

#print table
