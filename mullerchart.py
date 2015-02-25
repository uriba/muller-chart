import matplotlib.pyplot as mpl
import matplotlib.patches as pch
import numpy as np
from scipy.interpolate import spline
import pandas as pd
from pandas.io.parsers import read_csv
import pydot

freq_file_name = "Muller_Data.csv"      #The file name where the frequencies of the strains and timepoints are specified
hierarchy_filename = "Hierarchy.txt"    #The file that contains the strain hierarchy data in Python dictionary syntax with
                                        # keys being the strain names and values being a list of the daughter strains.

smoothing = False                   #should the resulting plot lines be smoothed via spline.


freq = read_csv(freq_file_name,header=None,index_col = 0)
freq = freq.transpose()
for col in freq.columns:
    freq.loc[freq[col]<0.005,col] = 0.0

# A dictionary stating for each strain its decendants. The order specified here determines the vertical order in which multiple decendants will be plotted. Earlier in the list = lower in the plot.
with open(hierarchy_filename,'r') as f:
    hierarchy = eval(f.read())

# plot the hierarchy tree:

graph = pydot.Dot(graph_type = 'graph')
def plot_tree(node,graphNode):
    for son in hierarchy[node]:
        sonNode = pydot.Node(son,label=son)
        graph.add_edge(pydot.Edge(graphNode,sonNode))
        plot_tree(son,sonNode)

root = pydot.Node("WT",label="WT")
plot_tree("WT",root)

graph.write_png("hierarchy.png")

# make sure each strain has abundance >= sum of decendent strains at every time point and print rounding performed.
def adjust_node(node,t):
    val = freq.loc[t,node]
    for son in hierarchy[node]:
        adjust_node(son,t)
    decendents_size = sum(freq.loc[t,hierarchy[node]])
    if decendents_size > val:
        print "adjusted %s at time %.1f by %.2f" % (node,t,decendents_size - val)
    newval = max(val,decendents_size)
    freq.loc[t,node] = newval
    return newval

for t in freq.index:
    adjust_node('WT',t)

# add time points for initiation of strains in cases where ancestor and decendent first arise in the same time point.
freq.set_index('Time',drop=False,inplace=True)

init_times = {'WT':0.0}

def initial_read_time(node):
    return min(freq.loc[freq[node] > 0,'Time'])    

def set_initiation_time(node):
    init_read_time = initial_read_time(node)
    if node not in init_times:
        init_time = max(freq.loc[freq['Time']<init_read_time,'Time'])
        init_times[node] = init_time
    init_time = init_times[node]
    for son in hierarchy[node]:
        if initial_read_time(son) == init_read_time:
            init_times[son] = (init_time + init_read_time)/2
            print "time added %f" % init_times[son]
        set_initiation_time(son)

set_initiation_time('WT')

newtimes = sorted(set(init_times.values()))
for i,t in enumerate(newtimes):
    if t not in freq['Time'].values:
        row = pd.DataFrame(index=[t])
        row['Time'] = t
        for node in hierarchy:
            if init_times[node] >= t:
                row[node] = 0.0
            else:
                start_time = max(freq.loc[freq['Time'] < t,'Time'])
                end_time = min(freq.loc[freq['Time'] > t,'Time'])
                start_val = freq.loc[start_time,node]
                end_val = freq.loc[end_time,node]
                time_frac = (t-start_time)/(end_time - start_time)
                row[node] = start_val+time_frac*(end_val - start_val)
        freq = freq.append(row)

times = sorted(freq.index)
freq.drop('Time',axis=1,inplace=True)
freq.sort_index(inplace=True)

pointbounds = {}  # pointbounds stores, for each strain and every time point, the beginning and ending vertical coordinates
                # of its slice.

#initialize the pointbounds dictionary
for i in hierarchy:
    pointbounds[i] = []

def calc_bounds(node,t,offset):
    size = freq.loc[t,node]
    pointbounds[node].append((offset,offset+size))
    sons = len(hierarchy[node])
    sons_size = sum(freq.loc[t,hierarchy[node]])
    interval = (size-sons_size)/(sons+1)
    offset += interval
    for i,son in enumerate(hierarchy[node]):
        calc_bounds(son,t,offset)
        offset+= interval + freq.loc[t,son]

for t in times:
    calc_bounds("WT",t,0) # For every time point recursively calculate the slices each strain occupies

# Assign colors to the different strains
colors = {"WT":["white","white"],
'0-xylE':[(0.8,0.8,1.0),(0.8,0.8,1.0)],
'0-topA':[(0.6,0.6,0.8),(0.6,0.6,0.8)],
'0-crp':[(0.4,0.4,0.6),(0.4,0.4,0.6)],
'0-yjiY':[(0.2,0.2,0.5),(0.2,0.2,0.5)],
'2-mlc+2':[(0.0,0.4,0.5),(0.0,0.0,0.0)],
'2-malE':[(0.0,0.6,0.6),(0.0,0.6,0.6)],
'2-thrA+2':[(0.0,0.75,0.5),(0.0,0.75,0.5)],
'2-prs+2':[(0.0,1.0,0.5),(0.0,1.0,0.5)],
'1-fliF':[(0.35,0.0,0.5),(0.0,0.0,0.0)],
'1-mlc+2':[(0.5,0.0,0.5),(0.5,0.0,0.5)],
'1-prs+7':[(0.5,0.0,0.75),(0.5,0.0,0.75)],
'1-cbdA':[(0.75,0.0,0.75),(0.75,0.0,0.75)], 
            'N-xylA*':["0.7","0.7"],
            'N-crp*':["0.5","0.5"],
            'N-rpoB*':[(0.0,0.5,0.3),(0.0,0.5,0.3)],
            'N-brnQ*':[(0.0,0.85,0.4),(0.0,0.85,0.4)], 
            'N-nadB*':[(0.0,0.4,0.3),(0.0,0.4,0.3)],
            'N-ptsI*':[(0.0,0.4,0.2),(0.0,0.4,0.2)],
            }
nodes = ["WT", '0-xylE','N-xylA*', 'N-crp*','0-topA', '0-crp', '0-yjiY', '1-fliF','2-mlc+2', '2-malE', '2-thrA+2', '2-prs+2', 'N-rpoB*',# "2-rpoB-malE",
 '1-mlc+2',
'1-prs+7', '1-cbdA','N-brnQ*','N-nadB*','N-ptsI*']

fig = mpl.figure(figsize = (14,6))
plt = fig.add_subplot(111)

# These are needed for the legend
handles = []
labels = []

# Loop through the strains and plot each one's slices.
def plot_node(node):
    coords = {'time':[],'ymin':[],'ymax':[]}
    tstart = None
    for t in range(len(times)):
        (ptmin,ptmax) = pointbounds[node][t]
            # For nicer visualization we omit slices of width 0 as they clutter the graph. We do need to include
            # slices of width zero if they either preceed or succeed a time point with that slice being non zero,
            # to show the emergence or decline of that strain.
        if ptmax-ptmin>0.0005:   #threshold to avoid rounding issues
            if tstart is None:
                tstart = times[max(0,t-1)]
            tend = times[min(t+1,len(times)-1)]
        coords["time"].append(times[t])
        coords["ymin"].append(ptmin)
        coords["ymax"].append(ptmax)
    if smoothing: #Smoothing can be applied to make the plot more visually appealing but with spline it does not always produce the desired results
        time = np.linspace(coords["time"][0],coords["time"][-1],100)
        ymin = spline(coords["time"],coords["ymin"],time)
        ymax = spline(coords["time"],coords["ymax"],time)
    else:
        time = np.array(coords['time'])
        ymin = np.array(coords['ymin'])
        ymax = np.array(coords['ymax'])
    indices = (time>=tstart) & (time<=tend)
    plt.fill_between(time[indices],ymin[indices],ymax[indices],facecolor=colors[node][0],edgecolor = colors[node][1],label=node)
    # Take care of the legend.
    handles.append(pch.Patch(facecolor = colors[node][0],edgecolor = "0.0",label = node))
    labels.append(node)
    #overlay decendents
    for son in hierarchy[node]:
        plot_node(son)

plot_node('WT')
plt.set_ylim(0,1.1)
plt.set_xlim(0,20.5)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='minor', labelsize=18)
plt.set_xlabel("time [weeks]", fontsize = 20)
plt.set_ylabel("population fraction", fontsize = 20)
mpl.figlegend(handles,labels,loc="right") 
mpl.subplots_adjust(right=0.8)
# And violla, our marvellous plot...        
mpl.savefig("muller-chart.pdf")
