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

randomSeed = int(time.time())
np.random.seed(randomSeed)
columns = 4
a_g = 0.0
b_g = 0.790391
c_g = -19.8928
a_r = 0.0
b_r = 0.766309
c_r = -19.0334
#SAMPLE f: sigma_g = a_g + np.exp(b_g*(y[i,2] + diff_mag) + c_g)

#d = raw_input('GC distance:')
#distance = eval(d)
#starfile = raw_input('Star file name:')
starfile = "5904_cluster_stars.txt"
outfile = 'convolved_5904_cluster_lilcon.txt'
distance = 7.5
far_d = 7.6
diff_mag = 5.0*(np.log10(far_d / distance))
print 'Convolving file', starfile, 'with random seed =', randomSeed

y = f.read_data(starfile)
l, w = y.shape
g = sc.zeros(l)
r = sc.zeros(l)
no_g = 0
no_r = 0
for i in range(l):
    #print y[i,2], y[i,3]
    sigma_g = (np.exp(b_g*(y[i,2]+diff_mag)) - np.exp(b_g*y[i,2]))*np.exp(c_g)
    sigma_r = (np.exp(b_r*(y[i,3]+diff_mag)) - np.exp(b_r*y[i,3]))*np.exp(c_r)
    if (sigma_g > 0.0):
        g[i] = np.random.normal(y[i,2], sigma_g)
    else:  
        g[i] = y[i,2]
        no_g = no_g + 1
    if (sigma_r > 0.0):
        r[i] = np.random.normal(y[i,3], sigma_r)
    else:  
        r[i] = y[i,3]
        no_r = no_r + 1
    #print sigma_g, sigma_r, g[i], r[i]
###Two columns, absolute g and g-r
print '-uncon g:', no_g, 'uncon r:', no_r
if columns == 2:
    m = sc.zeros((l,2))
    for i in range(l):
        m[i,0] = g[i] - 5.*(np.log10(distance*1000) - 1.)
        m[i,1] = g[i] - r[i]
###Four columns:  Ra, Dec, app_g, app_r; usually used for .csv output
elif columns == 4:
    m = sc.zeros((l,4))
    for i in range(l):
        m[i,0] = y[i,0]
        m[i,1] = y[i,1]
        m[i,2] = g[i]
        m[i,3] = r[i]
else:
    m = [0, 0, 0, 0]
    print "-Please use a valid column number"

if f.write_data(m, outfile, ' ') == 1:
    print '-Convolved data written successfully'


