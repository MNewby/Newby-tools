import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import files as f

'''convolves data sets - Matthew Newby '''

randomSeed = 151515
np.random.seed(randomSeed)
columns = 4
a_g = 0.0
b_g = 0.790391
c_g = -19.8928
a_r = 0.0
b_r = 0.766309
c_r = -19.0334

#d = raw_input('GC distance:')
#distance = eval(d)
#starfile = raw_input('Star file name:')
starfile = "5904_cluster_stars.txt"
distance = 7.5
far_d = 23.2
diff_mag = 5.0*(np.log10(far_d / distance))

y = f.read_data(starfile)
l, w = y.shape
g = sc.zeros(l)
r = sc.zeros(l)
for i in range(l):
    #print y[i,2], y[i,3]
    sigma_g = a_g + np.exp(b_g*(y[i,2] + diff_mag) + c_g)
    sigma_r = a_r + np.exp(b_r*(y[i,3] + diff_mag) + c_r) 
    g[i] = np.random.normal(y[i,2], sigma_g)
    r[i] = np.random.normal(y[i,3], sigma_r)
    #print sigma_g, sigma_r, g[i], r[i]
###Two columns, absolute g and g-r
if columns == 2:
    m = sc.zeros((l,2))
    for i in range(l):
        m[i,0] = g[i] - 5.*(np.log10(distance*1000) - 1.)
        m[i,1] = g[i] - r[i]
###Four columns:  Ra, Dec, app_g, app_r
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

if f.write_data(m, 'convolved_5904_cluster.csv', ',') == 1:
    print '-Convolved data written succesfully'


