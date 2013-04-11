#! /usr/bin/env python

import scipy as sc
import math as m
import numpy as np
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import gauss_fit_enc as ga

"""This program contains code for plotting globular cluster
  analyses.

Matthew Newby (RPI), June 5, 2010
"""

def plot_infiles(data_in, back_in, GCname='test', to_file=1):
    #plots data. data_in is 4-columns, ra,dec,g,r; 
    #to plot one file, pass it in twice.
    plt.figure()
    plt.scatter(data_in[:,0], data_in[:,1], c='white', marker='o')
    plt.scatter(back_in[:,0], back_in[:,1], c='black', marker='+')
    plt.axis('scaled')
    #plt.title('Sky Plot of cut data from ' + GCname)
    plt.xlim(np.ma.min(back_in[:,0]), np.ma.max(back_in[:,0]), fontsize='small')
    plt.ylim(np.ma.min(back_in[:,1]), np.ma.max(back_in[:,1]), fontsize='small')
    locs, lbls = plt.xticks()
    print locs, lbls
    labels = []
    for i in range(len(lbls)):  labels.append(str(locs[i]))
    plt.xticks(locs[1:-1], labels[1:-1], fontsize=8, rotation=0)
    plt.yticks(fontsize=8)
    plt.xlabel('RA', fontsize=8)
    plt.ylabel('Dec', fontsize=8)
    if to_file == 1:
        file_str = GCname + '_skyplot.ps'
        plt.savefig(file_str, papertype='letter')
        print '-skyplot created and saved as', file_str
    else:
        plt.show()   
        print '-skyplot created successfully'
    plt.close('all')
    return 1


def plot_hist(hist_in, GCname='test_hist', plot_fit=0, fit_params=[0.0,0.0,0.0,0.0],
    plot_params=0, to_file=1):
    #Make sure hist_in and params are in wanted format - abs or app magnitude.
    if plot_fit == 1:
        x_line = sc.arange( (np.ma.min(hist_in[:,0])*0.9), (np.ma.max(hist_in[:,0])*1.1), 0.01)
        y_line = sc.zeros(len(x_line))
        for i in range(len(x_line)):
            if x_line[i] < fit_params[0]:
                sigma = fit_params[1]
            else:
                sigma = fit_params[2]
            exponent = -0.5 * ( (x_line[i] - fit_params[0]) / sigma )**2
            stuff = m.exp(exponent)
            y_line[i] = stuff*fit_params[3]
        if plot_params == 1:
            txt_string = '$\mu_M=$ ' + str(fit_params[0][:6]) + '\n' \
            + '$\sigma_l=$ ' + str(fit_params[1])[:6] + '\n' \
            + '$\sigma_r=$ ' + str(fit_params[2])[:6] + '\n' \
            + '$A=$ ' + str(fit_params[3])[:6] + '\n' \
            + '$fwhm=$ '+ str(ga.full_width_half_max(fit_params[1],fit_params[2]))[:6]
        #Creates a line to simulate binning
    x_bins = sc.arange( (np.ma.min(hist_in[:,0])*0.9), (np.ma.max(hist_in[:,0])*1.1), 0.01)
    y_bins = sc.zeros(len(x_bins))
    span = hist_in[1,0] - hist_in[0,0]
    bin_min = hist_in[0,0] + (span/2.0) #Starts at first bin-bin changover
    j=0
    for i in range(len(x_bins)):
        if (x_bins[i] > ( bin_min + (span*j) )):
            j = j + 1
            if j > (len(hist_in[:,0])-1): j = (len(hist_in[:,0])-1)
        #print hist_in.shape, len(x_bins), len(y_bins)
        #print i, j
        y_bins[i] = hist_in[j,1]
    plt.figure(figsize=(8,6))
    if plot_fit == 1:
        plt.plot(x_line, y_line, 'g-')
        if plot_params == 1:
            plt.text(np.ma.min(x_line)+2.0, (np.ma.max(y_line)*0.6), (txt_string))
    plt.plot(x_bins, y_bins, 'b-')
    plt.scatter(hist_in[:,0], hist_in[:,1], marker='d' )
    #plt.title('Histogram of ' + GCname)
    plt.xlabel(r'$M_g$', fontsize=20)
    plt.ylabel('counts', fontsize=18)
    plt.xlim(0.0, 10.0, fontsize=14)
    plt.ylim((np.ma.min(y_bins)*0.5), (np.ma.max(y_bins)*1.2), fontsize=14)
    if to_file == 1:
        file_str = GCname + '_hist.ps'
        plt.savefig(file_str, papertype='letter')
        print '-histogram created and saved as', file_str
    else:
        plt.show()   
        print '-histogram created successfully'
    plt.close('all')
    return 1


def plot_HR_diagram(data_in, distance, GCname='test', ylimit=[1.0,7.0],
    xlimit=[0.0,0.6], to_file=1):
    #data_in should be 4-columns:  ra,dec,g,r
    length, width = data_in.shape
    gminusr = sc.zeros(length)
    abs_g = sc.zeros(length)
    for i in range(length):
        gminusr[i] = (data_in[i,2] - data_in[i,3])
        abs_g[i] = ( data_in[i,2] - 5.*(m.log10(distance*1000) - 1.) )
    plt.figure(figsize=(8,6))
    plt.scatter(gminusr, abs_g, marker='+')
    plt.title('HR diagram of ' + GCname)
    plt.xlabel('g-r')
    plt.ylabel('Mg')
    plt.xlim(xlimit[0], xlimit[1])
    plt.ylim(ylimit[1], ylimit[0])
    if to_file == 1:
        file_str = GCname + '_HR.ps'
        plt.savefig(file_str, papertype='letter')
        print '-HR diagram created and saved as', file_str
    else:
        plt.show()   
        print '-HR diagram created successfully'
    plt.close('all')
    return 1


def plot_HR_contour():
    return 0



