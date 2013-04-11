#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import files as f

"""This program plots data from the gc project

Matthew Newby (RPI), Jan 28, 2010
"""

#Data
data = f.read_data('best_results_firstpoint.txt')
names = ['NGC 4147', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 5466', 'NGC 5904', 'NGC 6205',
         'NGC 7078', 'NGC 7089', 'Pal 5'] #'NGC 6341', 
columns = ['Distance(kpc)', 'Fwhm', '$\mu$ (magnitudes)', '$\sigma_l$', '$\sigma_r$',
           'Amplitude (counts)', '$\mu$ error', '$\sigma_l$ error', '$\sigma_r$ error',
           'Amplitude error', 'Average Age(Gyr)', 'Fit Age(Gyr)', '<[Fe/H]> (Harris)',
           '<[Fe/H]> (CG97 Scale)', 'Harris Distance(kpc)']

#x values
i = 13
#y values
j = 3
#y_error values
e = j + 4

plt.figure(1)
#plot error bars!!!
plt.errorbar(data[:,i], data[:,j], yerr=data[:,e], fmt=None)
plt.scatter(data[:,i], data[:,j], marker='d')
#plt.ylim((np.ma.max(data[:,j])+0.5), (np.ma.min(data[:,j])-0.5))
for k in range(len(names)):
    plt.text(data[k,i]+0.05, data[k,j]+0.01, names[k], fontsize=10)
#plt.title('Plot of age vs metallicity, isocrone')
plt.xlabel(columns[i])
plt.ylabel(columns[j])
plt.show()
