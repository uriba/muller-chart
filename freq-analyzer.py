import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
from itertools import combinations

threshold = 0.05 #assumed measurement error.
freq_file_name = "groups.csv"      #Filename containing the frequencies of the strains and timepoints.

#Read the file
freq = read_csv(freq_file_name,header=None,index_col = 0)
freq = freq.transpose()
freq.index = freq["Time"]
del(freq["Time"])
freq = freq/100
print freq

#Iterate over all pairs of strains and determine their relation
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

