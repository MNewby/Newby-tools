import math as ma
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import astro_coordinates as co
import files as fi

rad = ma.pi/180.0
deg = 180.0/ma.pi

#data = fi.read_data("../../../Desktop/SgrTriax_DYN.dat")
data = fi.read_data("../../../Desktop/oblate.dat")
print "# - Data sucessfuly read in"
x, y = [], []
count = 0
p1 = [0.0, 160.0]
p2 = [160.0, -220.0]
p3 = [260.0, 60.0]
a1 = (p1[1]-p2[1])/(p1[0]-p2[0])
b1 = p1[1] - a1*p1[0]
a2 = (p2[1]-p3[1])/(p2[0]-p3[0])
b2 = p2[1] - a2*p2[0]
print a1, b1, a2, b2

for i in range(len(data[:,0])):
    #if data[i,0] > 150.0 and data[i,0] < 185.0:  count=count+1; continue  # lambda Cut
    #if abs(data[i,1]) > 1.0:  count=count+1; continue  # beta Cut
    """ For L05 models """
    if data[i,15] > 2.0:  count=count+1; continue
    #  Kill below first line
    if data[i,0] > p1[0] and data[i,0] < p2[0]:  
        if data[i,12] < (a1*data[i,0] + b1):  count=count+1; continue
    # kill above second line
    if data[i,0] > p2[0] and data[i,0] < p3[0]:  
        if data[i,12] > (a2*data[i,0] + b2):  count=count+1; continue
    #if data[i,13] > 72.4:  count=count+1; continue
    #if data[i,10] > 131.0:  count=count+1; continue
    #if data[i,12] < -105.0:  count=count+1; continue
    x.append(data[i,0])
    y.append(data[i,12])
    """ For L10 models 
    if data[i,24] > 4.0:  count=count+1; continue
    #if data[i,0] < 180.0 and data[i,25] > 0.0:  count=count+1; continue
    #if data[i,0] > 180.0 and data[i,25] < 0.0:  count=count+1; continue
    x.append(data[i,6])
    y.append(data[i,8])""" 
print count
plt.figure()
plt.scatter(x,y,s=1)
#plt.xlim(-100.0, 100.0)
#plt.ylim(-100.0, 100.0)
plt.show()
