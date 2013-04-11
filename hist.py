#! /usr/bin/env python  #'Bang' line - modify as needed

import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import sys

"""This script makes a 'responsible' histogram, outputs the bins and heights as
data, and plots the resulting histogram.

Matthew Newby (RPI), March 18, 2011
"""
#Maybe add drop-down lines from bins?

#Bins data into a histogram
def make_hist(data, bin_size, spread=[]):
    if spread == []:
        spread = [np.ma.min(data), np.ma.max(data)]
    nbins = int(m.ceil( (spread[1] - spread[0]) / bin_size ) )
    counts = sc.zeros((nbins,2), float)
    #Get bin centers
    for i in range(nbins):
        counts[i,0] = spread[0] + (i*bin_size) + (bin_size / 2.0)
    #add each point to the correct bin
    for i in range(len(data)):
        if (data[i] >= spread[1]): continue
        index = int( (data[i] - spread[0]) / bin_size)
        if (index < 0):  continue
        counts[index,1] = counts[index,1] + 1.0
    print '#-histogram created successfully'
    return counts

#Makes a histogram that includes the cumulative number of counts in each bin
def cumulative_hist(data, bin_size, spread=[]):
    if spread == []:
        spread = [np.ma.min(data), np.ma.max(data)]
    nbins = int(m.ceil( (spread[1] - spread[0]) / bin_size ) )
    counts = sc.zeros((nbins,2), float)
    #Get bin centers
    for i in range(nbins):
        counts[i,0] = spread[0] + (i*bin_size) + (bin_size / 2.0)
    #add each point to the correct bin
    for i in range(len(data)):
        if (data[i] > spread[1]): continue
        index = int( (data[i] - spread[0]) / bin_size)
        if (index < 0):  continue
        counts[index,1] = counts[index,1] + 1.0
    #Accumulate bins
    for i in range(1, nbins):
        counts[i,1] = counts[(i-1),1] + counts[i,1]
    print '#-Cumulative histogram created successfully;', sc.sum(counts[:,1]), 'items binned'
    return counts

#Takes in a 2-D array (bin centers and counts) and creates a plot.
def plot_histogram(hist_data, bin_size, limits=[], name='hist_out', x_label=''):
    edges = [(np.ma.min(hist_data[:,0])-(3.0*(bin_size/2.0)) ),
        (np.ma.max(hist_data[:,0])+(3.0*(bin_size/2.0)) )]
    if limits == []:
        limits = [(np.ma.min(hist_data[:,0])-(bin_size/2.0) ),
        (np.ma.max(hist_data[:,0])+(bin_size/2.0) )]
    nbins = len(hist_data[:,0])
    array_size = (nbins + 2)*2 # 2 points per bin, including extra spacer bins
    #Generate and populate array
    hist_array = sc.zeros((array_size,2), float)
    x, j, k = (edges[0]+bin_size), 0, 0
    #set x,y for two spots, move ahead, skip next turn
    for i in range(2, (array_size-2)):
        if (k == 0):
            hist_array[i,0] = x
            x = x + bin_size
            hist_array[(i+1),0] = x
            y = hist_data[j,1]
            hist_array[i,1] = y
            hist_array[(i+1),1] = y
            j = j + 1
            k = 1
        else:
            k = 0
            continue
    #First 2 points
    hist_array[0,0] = edges[0]
    hist_array[1,0] = edges[0] + bin_size
    #Last 2 points
    hist_array[-2,0] = x
    x = x + bin_size
    hist_array[-1,0] = x
    # Generate Plot
    plt.figure(1)
    plt.plot(hist_array[:,0], hist_array[:,1], 'b-')    
    plt.xlabel(x_label)
    plt.ylabel('counts', fontsize=10)
    plt.xlim(limits[0], limits[1])
    #save plot
    file_str = name + '_hist.ps'
    plt.savefig(file_str, papertype='letter')
    print '#-Histogram plotted and saved as', file_str
    plt.close('all')
    return 1

#plot multiple histograms on a common axis,
# bin_heights is a list of separate histogram series, of same length as common_x
def plot_multiple_hist(common_x, bin_heights, bin_size, limits=[], name='multi', x_label=''):
    edges = [(np.ma.min(common_x)-(3.0*(bin_size/2.0)) ),
        (np.ma.max(common_x)+(3.0*(bin_size/2.0)) )]
    if limits == []:
        limits = [(np.ma.min(common_x)-(bin_size/2.0) ),
        (np.ma.max(common_x)+(bin_size/2.0) )]
    nbins = len(common_x)
    # Build a list of values for each series (3 points each bin) to plot later
    series_list = []
    for i in range(len(bin_heights)):
        series = [0.0]  #first zero point
        for j in range(len(bin_heights[i])):
            series.append(bin_heights[i][j])
            series.append(bin_heights[i][j])
            series.append(0.0)
        series_list.append(series)
    # Build a matching list for the common axis
    x_axis = [common_x[0]-(bin_size/2.0)]  # first point
    for i in range(len(common_x)):
        x_axis.append(common_x[i]-(bin_size/2.0))
        x_axis.append(common_x[i]+(bin_size/2.0))
        x_axis.append(common_x[i]+(bin_size/2.0))
    # Plot this bad boy    
    plt.figure(1)
    for i in range(len(series_list)):
        plt.plot(x_axis, series_list[i])    
    plt.xlabel(x_label)
    plt.ylabel('counts', fontsize=10)
    plt.xlim(limits[0], limits[1])
    #save plot
    file_str = name + '_hist.ps'
    plt.savefig(file_str, papertype='letter')
    print '#-Multiple Histogram plotted and saved as', file_str
    plt.close('all')
    return 1

    
if __name__ == '__main__':
    print '#- number of arguments:', (len(sys.argv)-1)
    # Read arguments
    filename = sys.argv[1]
    column = int(sys.argv[2])
    size = float(sys.argv[3])
    if len(sys.argv) > 4:
        spread = [float(sys.argv[4]), float(sys.argv[5])]
    else:  spread = []
    if len(sys.argv) > 6:
        name = sys.argv[6]
    else:  name = 'quick'
    # Load Data
    suffix = filename[-4:]
    if suffix == '.txt':
        data = f.read_data(filename)
    elif suffix == '.csv':
        data = f.read_csv(filename)
    # Run scripts
    to_bin = data[:,column]
    reg_hist = make_hist(to_bin, size, spread)
    plot_histogram(reg_hist, size, name=(name+'_normal') )
    f.write_data(reg_hist, fileout=(name+'_normal.txt'), header='# Centers, Counts')
    cum_hist = cumulative_hist(to_bin, size, spread)
    plot_histogram(cum_hist, size, name=(name+'_cumulative') )
    f.write_data(cum_hist, fileout=(name+'_cumulative.txt'), header='# Centers, Counts')
    print '#done with quick histograms'