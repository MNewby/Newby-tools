import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import files as f

'''selects out only a portion of a data set
Used for f-turnoff stars right now.
---Matthew Newby '''

#Takes in a 4-column data file, with g and r as the last two columns,
# and returns only the stars in the given g-r range.
#could generalize to generic colors/column numbers.
def select_stars(data_in, low_limit = 0.0, high_limit = 0.3):
    length, width = data_in.shape
    data_index = sc.zeros(length, int)    
    for i in range(length):
        if (data_in[i,2] - data_in[i,3]) >= low_limit:
            if (data_in[i,2] - data_in[i,3]) <= high_limit:
                data_index[i] = 1
    s = sc.sum(data_index)
    print '-Stars in data set:', length, 'Stars in f-turnoff:', s
    new_data = sc.zeros((s,width))
    j = 0
    for i in range(length):
        if data_index[i] == 1:
            new_data[j,:] = data_in[i,:]
            j = j + 1
    if (j != s):  print '-!!! stars out not equal to f-turnoff:', j, s
    return new_data
