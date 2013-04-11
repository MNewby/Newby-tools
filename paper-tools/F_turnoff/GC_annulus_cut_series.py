#! /usr/bin/env python  #'Bang' line - modify as needed

import math as m
import numpy as np
import scipy as sc
import files as f
import gctools_enc as gce
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import sys

"""This is script takes a globular cluster, and cuts several annuli at differing radii.
The plots each annuli as a CMD and produces a histogram of g-magnitude values.
CMDs?  g vs. g-r AND u-g vs g-r?

Matthew Newby (RPI), Mar 15, 2011
"""

# Get Data
#Data is: ra, dec, g, r, u, flags
filename = 'noU_NGC_6205.csv'
name = '6205_test'
center_ra, center_dec = 250.423, 36.460 #NGC 6205
#center_ra, center_dec = 229.013, -0.123 #Pal 5
radii_step = 0.04
steps = 25

suffix = filename[-4:]
if suffix == '.txt':
    data = f.read_data(filename)
elif suffix == '.csv':
    data = f.read_csv(filename) 


'''Initial Skyplot'''
plt.figure(1)
data_1 = data[:,0]
data_2 = data[:,1]
#Strange indices are due to phantom liast point - not in data, but added by code somehow...
plt.scatter(data_1[:-1], data_2[:-1], 1, 'k', 'o')
plt.xlabel('ra')
plt.ylabel('dec')
plt.axis('scaled')
y_limits = plt.ylim()
x_limits = plt.xlim()
plot_file = name + '_initial_skyplot.ps'
plt.savefig(plot_file, papertype='letter')
plt.close('all')

#Loop over radii limits
for i in range(steps):
    inner_cut = i*radii_step
    outer_cut = (i+1)*radii_step
    annulus = gce.make_cut(data, center_ra, center_dec, inner_cut, outer_cut)
    if len(annulus[:,0]) == 0:
        print "#-No stars remained in this cut"
        continue
#plot CMDs and histograms
    '''Skyplot'''
    plt.figure(2)
    data_1 = annulus[:,0]
    data_2 = annulus[:,1]
    #Strange indices are due to phantom liast point - not in data, but added by code somehow...
    plt.scatter(data_1[:-1], data_2[:-1], 1, 'k', 'o')
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.axis('scaled')
    plt.xlim(x_limits[0], x_limits[1])
    plt.ylim(y_limits[0], y_limits[1])
    plot_file = name + '_' + str(i) + '_skyplot.ps'
    plt.savefig(plot_file, papertype='letter')
    plt.close('all')
    '''Plot g vs. g-r'''
    plt.figure(3)
    data_g = annulus[:,2]
    data_r = annulus[:,3]
    new_x = (data_g - data_r)
    plt.scatter(new_x, data_g, 1, 'k', 'o')
    #y_limits = plt.ylim()
    #plt.ylim(limits[1], limits[0])
    plt.ylim(25.0, 10.0)
    plt.xlim(-0.3, 0.6)
    plt.xlabel(r'$(g-r)_0$')
    plt.ylabel(r'$g_0$')
    plot_file = name + '_' + str(i) + '_cmd.ps'
    plt.savefig(plot_file, papertype='letter')
    plt.close('all')
    '''Plot u-g vs g-r'''
    plt.figure(4)
    data_g = annulus[:,2]
    data_r = annulus[:,3]
    data_u = annulus[:,4]
    new_x = (data_g - data_r)
    new_y = (data_u - data_g)
    plt.scatter(new_x, new_y, 1, 'k', 'o')
    plt.xlabel(r'$(g-r)_0$')
    plt.ylabel(r'$(u-g)_0$')
    plt.ylim(0.0, 4.0)
    plot_file = name + '_' + str(i) + '_cc.ps'
    plt.savefig(plot_file, papertype='letter')
    plt.close('all')
    '''g histogram'''
    plt.figure(5)
    data_g = annulus[:,2]
    nbins = 20
    plt.hist(data_g, nbins)  #20
    plt.xlabel(r'$g_0$')
    #plt.ylabel('counts')
    plt.xlim(10,30)
    plot_file = name + '_' + str(i) + '_hist.ps'
    plt.savefig(plot_file, papertype='letter')
    plt.close('all')
    print "#-Finished with iteration", i
print '#-Done'