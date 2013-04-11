import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import files as f
import time

'''convolves GC data sets - 
---Starfile should be the full HR color-magnitude diagram stars; as the
   convolution will push some stars into/out of the f-turnoff bin.
   Format should be: 'RA(space)Dec(space)dered_g(space)dered_r(space)'
   With one star per line.
---Maximum distance at 23.2
---Steps of 2.58 for 7078
---Steps of 1.9 for 6205
---Steps of 1.95 for 7098
---ERROR somewhere!  using same distances does not produce and unconvolved data set!
Matthew Newby, RPI, 2010 '''

#Takes a 4 column data set in, with g in the 3rd column and r in the 4th column.
#Convolves r and g to the SDSS errors at a new distance, 'con_dist'
#Returns a 4 column data set, with columns 1 & 2 unchanged, and g&r convolved.
#Possibly make this more generic, so g and r can be any column? 
def convolve(data_in, real_dist, con_dist):
    randomSeed = int(time.time())
    np.random.seed(randomSeed)
    a_g = 0.0
    b_g = 0.790391
    c_g = -19.8928
    a_r = 0.0
    b_r = 0.766309
    c_r = -19.0334
    diff_mag = 5.0*(np.log10(con_dist / real_dist))
    print 'Convolving file', starfile, 'with random seed =', randomSeed
    length, width = data_in.shape
    g = sc.zeros(length)
    r = sc.zeros(length)
    for i in range(length):
        sigma_g = (np.exp(b_g*(data_in[i,2]+diff_mag)) - np.exp(b_g*data_in[i,2]))*np.exp(c_g)
        sigma_r = (np.exp(b_r*(data_in[i,3]+diff_mag)) - np.exp(b_r*data_in[i,3]))*np.exp(c_r)
        if (sigma_g > 0.0):
            g[i] = np.random.normal(data_in[i,2], sigma_g)
        else:  
            g[i] = data_in[i,2]
        if (sigma_r > 0.0):
            r[i] = np.random.normal(data_in[i,3], sigma_r)
        else:  
            r[i] = data_in[i,3]
    m = sc.zeros((length,4))
    for i in range(length):
        m[i,0] = data_in[i,0]
        m[i,1] = data_in[i,1]
        m[i,2] = g[i]
        m[i,3] = r[i]
    return m