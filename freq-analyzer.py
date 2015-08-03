import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
from itertools import combinations

threshold = 0.05

freq_file_name = "Muller_Data.csv"      #The file name where the frequencies of the strains and timepoints are specified

freq = read_csv(freq_file_name,header=None,index_col = 0)
freq = freq.transpose()
for col in freq.columns:
    for i,row in enumerate(freq.index):
        prev = 0 if i == 0 else freq.iloc[i-1][col]
        following = 0 if i == len(freq.index)-1 else freq.iloc[i+1][col]
        if 0 < freq.loc[row,col] < 0.005 and prev == 0 and following == 0:
            print "rounded %s at time %d by %f" % (col,row,freq.loc[row,col])
            freq.loc[row,col] = 0.0
freq.index = freq["Time"]
del(freq["Time"])
print freq

for (x,y) in combinations(freq.columns,2):
    maybeSiblings = not max(freq[x]+freq[y])>1+threshold
    if((freq[x] < freq[y] + threshold).all()):
        print("%s %s %s" % (x,"<=" if maybeSiblings else "<",y))
    elif((freq[y] < freq[x] + threshold).all()):
        print("%s %s %s" % (y,"<=" if maybeSiblings else "<",x))
    elif(maybeSiblings):
        print("%s = %s" % (x,y))
    else:
        print("Inconsistent %s,%s" % (x,y))

