#! /usr/bin/env python

import files as f
import scipy as sc
import math as m
import numpy as np
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt


"""This program plots GC fit results as a ganged plot

Matthew Newby (RPI), Dec 11, 2010
"""

save = 1
save_name = 'GC_results_combined_metal.ps'
data = f.read_data('results_new_corrected_singles.txt')
names = ['NGC 4147', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 5466', 'NGC 5904', 'NGC 6205',
         'NGC 7078', 'NGC 7089', 'Pal 5'] #'NGC 6341',
x_column = 13  #13=CG97 metallicity, 11=fit age, 0 = distance
# !!! Order for Offsets:   sig_r, sig_l, mu, turnoff
""" Metallicity Offsets """
x_offset = [ [0.01, 0.01, 0.01, 0.01],
    [0.01, 0.01, 0.01, 0.01],
    [0.01, 0.01, 0.01, 0.01],
    [-0.11, -0.1, -0.1, -0.1],
    [-0.1, 0.01, 0.01, -0.1],
    [0.01, 0.01, 0.01, 0.01],
    [-0.02, -0.1, 0.01, 0.01],
    [0.01, 0.01, 0.01, 0.01],
    [0.01, 0.01, 0.01, 0.01],
    [0.01, 0.01, 0.01, 0.01] ]
y_offset = [ [0.0, 0.0, 0.0, -0.02],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.01],
    [-0.06, 0.0, -0.1, 0.0],
    [0.0, 0.0, 0.0, -0.03],
    [0.0, 0.0, 0.0, 0.0],
    [0.13, -0.05, -0.2, -0.08],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0] ]
"""  Age Offsets 
x_offset = [ [0.05, 0.05, 0.05, 0.05],
    [0.05, 0.05, 0.05, 0.05],
    [0.05, 0.05, 0.05, 0.05],
    [0.05, 0.05, 0.05, 0.05],
    [0.05, 0.05, 0.05, 0.05],
    [0.05, 0.05, 0.05, 0.05],
    [0.05, 0.05, 0.05, 0.05],
    [0.05, -0.4, -0.4, -0.4],
    [-0.35, -0.37, 0.05, 0.05],
    [0.05, 0.05, 0.05, 0.05] ]
y_offset = [ [0.0, 0.0, 0.0, 0.03],
    [0.0, 0.0, 0.0, 0.0],
    [-0.1, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, -0.1, 0.02],
    [0.0, 0.0, 0.0, 0.0],
    [-0.1, -0.05, 0.0, 0.04],
    [0.0, 0.0, 0.0, 0.0] ]"""
"""  Distance Offsets """

#ghost_point = sc.array( [[4.67,0.0736,-2.30],[0.753,0.0608,-2.30],[0.985,0.071,-2.30]] )


"""Plot Initializations"""
fig = plt.figure()
plt.subplots_adjust(hspace=0.001)

"""Iterate over parameters"""
subs = []
for i in [15,4,3,2]:
    if (i != 15):
        sp = plt.subplot(4,1,(i-1), sharex=subs[0])
    else:
        sp = plt.subplot(4,1,4)
    subs.append(sp)
    
    """Plot data and text"""
    if (i != 15):
        plt.errorbar(data[:,x_column], data[:,i], yerr=data[:,(i+4)], fmt='s',
                 color='black', ecolor='k', ms=3)
    else:
        plt.errorbar(data[:,x_column], data[:,i], yerr=data[:,17], fmt='s',
                 color='black', ecolor='k', ms=3)
    #For ghost point
    #plt.errorbar(ghost_point[i-2][2], ghost_point[i-2][0], yerr=ghost_point[i-2][1], fmt='o',
    #             color='white', ecolor='k', ms=4)
    #plt.text( (ghost_point[i-2][2]+0.02), (ghost_point[i-2][0]),'NGC 5053, C2009', fontsize=4)
    for j in range(len(names)):
        if (i != 15):
            plt.text( (data[j,x_column]+x_offset[j][abs(i-4)]), (data[j,i]+y_offset[j][abs(i-4)]),
                names[j], fontsize=6)
        else:
            plt.text( (data[j,x_column]+x_offset[j][3]), (data[j,i]+y_offset[j][3]),
                names[j], fontsize=6)


    """Setup axis labels and ticks"""
    if (x_column == 0):
        pass
        #plt.xlim(9.0, 14.0)  #Distance
        #plt.xticks(sc.arange(9.0, 14.0, 0.5), fontsize=8)  #Distance
    if (x_column == 11):
        plt.xlim(9.0, 14.0)  #Age
        plt.xticks(sc.arange(9.0, 14.0, 0.5), fontsize=8)  #Age
    if (x_column == 13):
        plt.xlim(-2.4, -1.0)  #Metallicity
        plt.xticks(sc.arange(-2.4, -0.8, 0.2), fontsize=8) #Metallicity
    if (i != 15):
        plt.setp(sp.get_xticklabels(), visible=False)
    else:
        if (x_column == 0):  plt.xlabel('Distance, kpc', fontsize=10)  #Distance
        if (x_column == 11):  plt.xlabel('Age, Gyr', fontsize=10)  #Age
        if (x_column == 13):  plt.xlabel('$<[Fe/H]>$ (CG97 scale)', fontsize=10) #Metallicity
    if (i == 4):
        plt.ylabel('$\sigma_r$ (magnitudes)', fontsize=8)
        plt.yticks(sc.arange(0.6, 1.5, 0.2), fontsize=8)
    elif (i == 3):
        plt.ylabel('$\sigma_l$ (magnitudes)', fontsize=8)
        plt.yticks(sc.arange(0.2, 1.1, 0.2), fontsize=8)
    elif (i == 2):
        plt.ylabel('$\mu$ (magnitudes)', fontsize=8)
        plt.yticks(sc.arange(4.0, 5.0, 0.2), fontsize=8)
    elif (i== 15):
        plt.ylabel('$(g-r)_0$ turnoff', fontsize=8)
        plt.yticks(sc.arange(0.0, 0.5, 0.1), fontsize=8)
        
"""Output"""
#if save == 1:
#    plt.savefig(save_name, papertype='letter')
#    print '#---data plot saved as', save_name
#else:
plt.show()
plt.close('all')
print '#---Done with data plot'