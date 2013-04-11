#! /usr/bin/env python

import files as f
import hist as h
import scipy as sc
import math as m
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


"""This script is for simple scripts.

Matthew Newby (RPI), May 4, 2011
"""

spread = 0.5

datafile = ['noU_NGC_4147_cluster.csv',
            'noU_NGC_5024_cluster.csv',
            'noU_NGC_5053_cluster.csv',
            'noU_NGC_5272_cluster.csv',
            'noU_NGC_5466_cluster.csv',
            'noU_NGC_5904_cluster.csv',
            'noU_NGC_6205_cluster.csv',
            'noU_NGC_7078_cluster.csv',
            'noU_NGC_7089_cluster.csv',
            'noU_Pal5_cluster.csv',
           ]
DISTANCE = [19.3, 18.7, 18.5, 10.4, 15.6, 8.0, 7.7, 11.0, 11.5, 21.0]

results_file = 'results_new_corrected_singles.txt'

results = f.read_data(results_file)

for i in range(len(datafile)):
    data = f.read_data(datafile[i], ',')
    mu = results[i,2] + 5.*(m.log10(DISTANCE[i]*1000) - 1.)
    holder=[]
    for j in range(len(data[:,0])):
        if data[j,2] > (mu+spread):  continue
        if data[j,2] < (mu-spread):  continue
        holder.append(data[j,2]-data[j,3])
    histogram = h.make_hist(holder, 0.01)
    h.plot_histogram(histogram, 0.01, name=('hist_'+str(i)), x_label=r'g-r')