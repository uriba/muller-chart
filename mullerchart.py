import matplotlib.pyplot as mpl
import matplotlib.patches as pch
import numpy as np
from scipy.interpolate import spline
import pandas as pd

normalize = True
smoothing = False
#the abundances of every strain as a function of time.
mutfreq_orig = {
"WT":[1,1,1,1,1,1,1,1,1,1,1,1],
"0-crp":[0,0.118881119,0.599878429,0.729227941,0.913538179,0.845454546,0.898015356,0.871893086,0.881793479,0.869071573,0.924837071,0.933035715],
"0-topA":[0,0.944954128,0.987012987,1,0.956989247,0.983050847,0.941176471,0.951219512,0.902654867,0.808988764,1,1],
"0-xylE":[0.33974359,0.971014493,1,1,0.988023952,0.950980392,0.926829268,0.973856209,0.955555556,0.886792453,0.979865772,0.985507246],
"0-yjiY":[0,0.020698443,0.053251186,0.150121938,0.965277778,0.958000887,0.935779817,0.937282135,0.936507937,0.865771812,0.985915493,0.991150442],
"1-malE":[0,0.008264463,0,0,0.114754098,0.234042553,0.369747899,0.315068493,0.820754717,0.735042735,0.114503817,0.380165289],
"1-mlc+2":[0,0,0,0.002666667,0.84365705,0.7241301,0.747863248,0.916050042,0.898124854,0.755029262,0.211945371,0.336181978],
"1-prs+2":[0.001980198,0,0.001574803,0,0,0,0,0,0.028377293,0.009381545,0.144418747,0.345030256],
"1-thrA+2":[0,0,0,0,0.01,0.17,0.02,0.01,0.34358872,0.726033412,0.24,0.39],
"2-cbdA":[0,0,0,0,0,0,0,0,0,0.01010101,0.113043478,0.044247788],
"2-fli":[0,0,0.011235955,0.125,0.024590164,0.060869565,0,0.024539877,0.007352941,0.234848485,0.789473684,0.63559322],
"2-icd":[0,0,0.061327832,0.06027668,0.111763237,0.049136391,0.359432347,0.294515064,0.383649456,0.207464147,0.670754717,0.442550505],
"2-mlc+2":[0,0,0,0,0,0.001494768,0,0,0.006483338,0.114128627,0.68983125,0.487441085],
"2-prs+7":[0,0,0,0.003430938,0,0.001347709,0,0,0.001335113,0.104886519,0.709374974,0.527642027],
"N-crp*":[0,0.057569745,0.517433467,0.56119193,0.00204499,0,0.005509642,0,0.004694836,0.004385965,0.002252252,0.002873563],
"N-xylA*":[0.047058824,0.66,0.052419355,0.012722646,0,0,0,0,0.002197802,0,0,0]
}


times = [2,4,5,6,8,10,11,13,16,18,19,20.5]
mutfreq = {}
for key in mutfreq_orig:
    mutfreq[key] = np.around(mutfreq_orig[key],2)
df = pd.DataFrame
df = df.from_dict(mutfreq,orient="index")
df = np.round(df,2)
df.columns = times

print df

#This dictionary states for each strain its decendants.
hierarchy = {"WT":['0-xylE'],
            '0-xylE':['0-topA'],
            '0-topA':['N-xylA*','0-crp'],
            '0-crp':['0-yjiY','N-crp*'],
            '0-yjiY':['1-mlc+2','2-icd'],
            '1-mlc+2':['1-malE'],
            '1-malE':['1-thrA+2'],
            '1-thrA+2':['1-prs+2'],
            '1-prs+2':[],
            '2-icd':['2-fli'],
            '2-fli':['2-mlc+2'],
            '2-mlc+2':['2-prs+7'],
            '2-prs+7':['2-cbdA'],
            '2-cbdA':[],
            'N-xylA*':[],
            'N-crp*':[],
            }

#the times used - in this example I simply used the index of the measurements.
print times

abundances = {}
for node in hierarchy:
    abundances[node]=[]

def set_abundances(node,t):
    sz = 0
    for son in hierarchy[node]:
        son_size = set_abundances(son,t)
        sz += son_size
    my_size = max(mutfreq[node][t],sz)
    if sz > mutfreq[node][t]:
        print "%s, at time %.1f, rounded %.2f" % (node,times[t],sz-mutfreq[node][t])
    abundances[node].append(my_size-sz)
    return my_size

for t in range(len(times)):            
    set_abundances("WT",t)
    if normalize:
        tot_size = 0
        for node in hierarchy:
            tot_size += abundances[node][t]
        scale = 1.0/tot_size
        for node in hierarchy:
            abundances[node][t] = abundances[node][t] * scale

    
for node in abundances:
    print node
    print abundances[node]

sizes = {} #sizes stores, for each strain, the size each "slice" of it occupies at every time point.
pointabdc = {}  #pointabdc stores, for each strain and every time point, the beginning and ending vertical coordinates
                # of each of its slices.

#initialize the sizes dictionary
for i in hierarchy:
    dec = len(hierarchy[i])
    ser = []
    for j in abundances[i]:
        ser.append(float(j)/(dec+1))
    sizes[i]=ser
    pointabdc[i] = []

#two recursive functions needed to initialize the pointabdc dictionary
def calc_size(node,t):
    nodesize = abundances[node][t]
    for son in hierarchy[node]:
        nodesize+=calc_size(son,t)
    return nodesize

def calc_splits(node,t,offset):
    slice_size = sizes[node][t] 
    points = [offset,offset+slice_size]
    for son in hierarchy[node]:
        calc_splits(son,t,points[-1])
        sz = calc_size(son,t)
        points.append(points[-1]+sz)
        points.append(points[-1]+slice_size)
    pointabdc[node].append(points)
    
for t in range(len(times)):
    calc_splits("WT",t,0) #for every time point do a recursive calculation of the slices each strain occupies

### generating the plotting data    
colors = {"WT":"0.3",
            '0-xylE':"#500000",
            '0-topA':"#700000",
            '0-crp':"#900000",
            '0-yjiY':"#a00000",
            '1-mlc+2':"#005000",
            '1-malE':"#007000",
            '1-thrA+2':"#009000",
            '1-prs+2':"#00c000",
            '2-icd':"#000030",
            '2-fli':"#000060",
            '2-mlc+2':"#000090",
            '2-prs+7':"#008fb2",
            '2-cbdA':"#00ccff",
            'N-xylA*':"0.5",
            'N-crp*':"0.7"
            }

fig = mpl.figure(figsize = (12,6))
plt = fig.add_subplot(111)
handles = []
labels = []
for node in sorted(pointabdc.keys()):
    for i in range(len(pointabdc[node][0])/2):
        coords = {'time':[],'ymin':[],'ymax':[]}
        for t in range(len(times)):
            ptmin = pointabdc[node][t][2*i]
            ptmax = pointabdc[node][t][2*i+1]
            if t == 0:
                diffprev = 0
            else:
                diffprev = pointabdc[node][t-1][2*i+1] - pointabdc[node][t-1][2*i]
            if t == len(times)-1:
                diffnext = 0
            else:
                diffnext = pointabdc[node][t+1][2*i+1] - pointabdc[node][t+1][2*i]
            if ptmax-ptmin+diffprev+diffnext > 0.005:
                coords["time"].append(times[t])
                coords["ymin"].append(ptmin)
                coords["ymax"].append(ptmax)
        if smoothing:
            softtimes = np.linspace(coords["time"][0],coords["time"][-1],100)
            softymin = spline(coords["time"],coords["ymin"],softtimes,order=3)
            softymax = spline(coords["time"],coords["ymax"],softtimes,order=3)
            plt.fill_between(softtimes,softymin,softymax,color=colors[node])
        else:
            plt.fill_between(coords['time'],coords['ymin'],coords['ymax'],color=colors[node],label=node)
    handles.append(pch.Patch(color = colors[node],label = node))
    labels.append(node)
plt.set_ylim(0,1.5)
plt.set_xlim(2,20.5)
mpl.figlegend(handles,labels,loc="upper right") 
mpl.subplots_adjust(right=0.8)
        
mpl.savefig("muller-demo2.pdf")
