import matplotlib.pyplot as mpl
import matplotlib.patches as pch
import numpy as np
from scipy.interpolate import spline
import pandas as pd
from pandas.io.parsers import read_csv

freq_file_name = "Muller_Data.csv"      #The file name where the frequencies of the strains and timepoints are specified
hierarchy_filename = "Hierarchy.txt"    #The file that contains the strain hierarchy data in Python dictionary syntax with
                                        # keys being the strain names and values being a list of the daughter strains.

normalize = False                   #should the total abundances be normalized to unity.
smoothing = False                   #should the resulting plot lines be smoothed via spline.

freq = read_csv(freq_file_name,header=None,index_col = 0)
freq = freq.transpose()
times = freq['Time']
freq = freq.drop('Time',1)
for col in freq.columns:
    freq[col] = np.around(freq[col],2)
freq['0-topA'] = freq['0-topA']-freq['N-xylA*'] - freq['N-crp*']
freq['0-crp'] = freq['0-crp']-freq['N-crp*']
freq['0-topA-tmp'] = 0.0
freq['fake1'] = 0.0
freq['fake2'] = 0.0

# A dictionary stating for each strain its decendants. The order specified here determines the vertical order in which multiple decendants will be plotted. Earlier in the list = lower in the plot.
hierarchy = eval(open(hierarchy_filename,'r').read())

# The abundances dictionary states the abundance of every strain in the population at every time plot. It is calculated
# from the mutfreq dictionary above by subtracting from the frequency of every mutation the frequencies of all of its
# sub-mutations (recursively), at each time point.
abundances = {}
for node in hierarchy:
    abundances[node]=[]

def set_abundances(node,t):
    sz = 0
    for son in hierarchy[node]:
        son_size = set_abundances(son,t)
        sz += son_size
    my_size = max(freq.loc[t,node],sz)
    if sz > freq.loc[t,node]: # If the data is contradictory, having a mutation the decendents of which exceeding its own frequency, print a report about it and the amount that it was rounded up by.
        print "%s, at time %.1f, rounded %.2f" % (node,times[t],sz-freq.loc[t,node])
    abundances[node].append(my_size-sz)
    return my_size

for t in range(1,len(times)+1):            
    set_abundances("WT",t)
    if normalize: # normalize the sum of abundances of all the strains to unity, which is not always the case in real experimental data.
        tot_size = 0
        for node in hierarchy:
            tot_size += abundances[node][t]
        scale = 1.0/tot_size
        for node in hierarchy:
            abundances[node][t] = abundances[node][t] * scale

    
sizes = {}  # Sizes stores, for each strain, the size each "slice" of it occupies at every time point (slice is a vertical
            # portion of the graph that is colored in that strain's color. At a given time point a strain may have few 
            # slices as they "wrap" every one of its decendent strains).

pointabdc = {}  # pointabdc stores, for each strain and every time point, the beginning and ending vertical coordinates
                # of each of its slices.

#initialize the sizes dictionary
for i in hierarchy:
    decendents = len(hierarchy[i])
    sizes[i] = []
    for j in abundances[i]:
        sizes[i].append(float(j)/(decendents+1))    # The width of every slice is the abundance of the strain divided 
                                                    # by its number of decendents + 1 so enough slices exist to wrap all the decendents.
    pointabdc[i] = []

# Two recursive functions needed to initialize the pointabdc dictionary:
# calc_size calculates, given a strain and a time point, the fraction of the population that node and its decendents occupy (recursively).
def calc_size(node,t):
    nodesize = abundances[node][t]
    for son in hierarchy[node]:
        nodesize+=calc_size(son,t)
    return nodesize

# calc_splits calculates the vertical beginning and ending point of each slice by interleaving them with the daughter strains of the given node (recursively).
def calc_splits(node,t,offset):
    slice_size = sizes[node][t] 
    points = [offset,offset+slice_size] # first slice starts at the starting offset of the node and is slice_size tall (as are all the slices)
    for son in hierarchy[node]:
        calc_splits(son,t,points[-1])
        sz = calc_size(son,t)
        points.append(points[-1]+sz)    # the next slice starts after the daughter strain ends
        points.append(points[-1]+slice_size) # and is again slice_size wide
    pointabdc[node].append(points)
    
for t in range(len(times)):
    calc_splits("WT",t,0) # For every time point recursively calculate the slices each strain occupies

# Assign colors to the different strains
colors = {"WT":"white",
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
            'N-crp*':"0.7",
            'N-rpoB*':"0.0",
            "1-rpoB-malE":"0.1",
            '0-topA-tmp':"0.0",
            'fake2':"0.0",
            'fake1':"0.0"
            }

fig = mpl.figure(figsize = (12,6))
plt = fig.add_subplot(111)

# These are needed for the legend
handles = []
labels = []

# Loop through the strains and plot each one's slices.
for node in sorted(pointabdc.keys()):
    for i in range(len(pointabdc[node][0])/2): # Each slice is defined by two lines - the lower and upper bounds of the slice.
        coords = {'time':[],'ymin':[],'ymax':[]}
        tstart = None
        for t in range(len(times)):
            ptmin = pointabdc[node][t][2*i]
            ptmax = pointabdc[node][t][2*i+1]
                # For nicer visualization we omit slices of width 0 as they clutter the graph. We do need to include
                # slices of width zero if they either preceed or succeed a time point with that slice being non zero,
                # to show the emergence or decline of that strain.
            if ptmax-ptmin>0.005:   #threshold to avoid rounding issues
                if not tstart:
                    tstart = times[max(1,t)]
                tend = times[min(t+2,len(times))]
            coords["time"].append(times[t+1])
            coords["ymin"].append(ptmin)
            coords["ymax"].append(ptmax)
        if not tstart:
            continue
        if smoothing: #Smoothing can be applied to make the plot more visually appealing but with spline it does not always produce the desired results
            time = np.linspace(coords["time"][0],coords["time"][-1],100)
            ymin = spline(coords["time"],coords["ymin"],time)
            ymax = spline(coords["time"],coords["ymax"],time)
        else:
            time = np.array(coords['time'])
            ymin = np.array(coords['ymin'])
            ymax = np.array(coords['ymax'])
        indices = (time>=tstart) & (time<=tend)
        plt.fill_between(time[indices],ymin[indices],ymax[indices],color=colors[node],label=node)
    # Take care of the legend.
    handles.append(pch.Patch(color = colors[node],label = node))
    labels.append(node)
plt.set_ylim(0,1.5)
plt.set_xlim(2,20.5)
mpl.figlegend(handles,labels,loc="upper right") 
mpl.subplots_adjust(right=0.8)
# And violla, our marvellous plot...        
mpl.savefig("muller-chart.pdf")
