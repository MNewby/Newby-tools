import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import files as f

'''selects out only a portion of a data set
Used for f-turnoff stars right now.
---Matthew Newby '''

starfile = "convolved_5904_cluster_lilcon.txt"
outfile = 'convolved_5904_cluster_lilcon_f.txt'
low_limit = 0.0
high_limit = 0.3

y = f.read_data(starfile)
l, w = y.shape
index = sc.zeros(l, int)
for i in range(l):
    if (y[i,2] - y[i,3]) >= low_limit:
        if (y[i,2] - y[i,3]) <= high_limit:
            index[i] = 1
s = sc.sum(index)
print '-Stars in data set:', l, 'Stars in f-turnoff:', s
new_data = sc.zeros((s,w))
j = 0
for i in range(l):
    if index[i] == 1:
        new_data[j,:] = y[i,:]
        j = j + 1
if (j != s):  print '-!!! stars out not equal to f-turnoff:', j, s
if f.write_data(new_data, outfile, ' ') == 1:
    print '-F-turnoff stars written successfully'
