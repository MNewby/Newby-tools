#! /usr/bin/env python

import sys
import files as f
import scipy as sc
import math as m
import numpy as np
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt


"""Plots the results from convolve analysis script "GC_CAscript.py"

Matthew Newby (RPI), Feb 16, 2011
"""

#CA_6205_1kpc_1kpc_100avg.txt
def plot_stuff(args):
    g_avg = 4.18
    filename=args[0]
    data = f.read_data(filename)
    l, w = data.shape
    blue_base, yellow_base, red_base = data[0,1], data[0,5], data[0,9]
    """ 2nd file stuff """
    """file2 = args[1]
    data2 = f.read_data(file2)
    l2, w2 = data2.shape
    blue_base2, yellow_base2, red_base2 = data2[0,1], data2[0,5], data2[0,9]"""
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    '''Plot Blue Box'''    
#    plt.subplot(3,1,1)
#    plt.plot( data[:,0], (data[:,1]/blue_base), 'b-')
#    plt.plot( data[:,0], (data[:,2]/blue_base), 'y-')
#    plt.plot( data[:,0], (data[:,3]/blue_base), 'r-')
    '''Plot Yellow Box'''    
#    plt.subplot(3,1,2)
#    plt.plot( data[:,0], (data[:,4]/yellow_base), 'b-')
#    plt.plot( data[:,0], (data[:,5]/yellow_base), 'y-')
#    plt.plot( data[:,0], (data[:,6]/yellow_base), 'r-')
    '''Plot Red Box'''    
#    plt.subplot(3,1,3)
#    plt.plot( data[:,0], (data[:,7]/red_base), 'b-')
#    plt.plot( data[:,0], (data[:,8]/red_base), 'y-')
#    plt.plot( data[:,0], (data[:,9]/red_base), 'r-')
    '''Plot Only turnoff stars'''
    ax1.plot( data[:,0], (data[:,5]/yellow_base), 'y-')
    ax1.plot( data[:,0], (data[:,2]/yellow_base), 'b-')
    ax1.plot( data[:,0], (data[:,8]/yellow_base), 'r-')
    ax1.plot( data[:,0], ((data[:,5]+data[:,2]+data[:,8])/yellow_base), 'k-')
    ''' Agiain!  For files '''
    """ax1.plot( data2[:,0], (data2[:,5]/yellow_base), 'y:')
    ax1.plot( data2[:,0], (data2[:,2]/yellow_base), 'b:')
    ax1.plot( data2[:,0], (data2[:,8]/yellow_base), 'r:')
    ax1.plot( data2[:,0], ((data2[:,5]+data2[:,2]+data2[:,8])/yellow_base), 'k:')"""
    '''Fancy double x-axis plotting stuff'''
    x_tick = ax1.get_xticks()
    ax1.set_xticks(x_tick)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel(r'$d_{\rm eff}$, kpc', fontsize=10)
    plt.ylabel(r'Fraction of Initial', fontsize=10)
    new_ticks = []
    for i in range(len(x_tick)):
        if x_tick[i] == 0:
            new_ticks.append("") #(g_avg)
        else:
            new_ticks.append(round( (5.0*m.log10(x_tick[i]/0.01) + g_avg), 1) )
    x_lim = ax1.get_xlim()
    new_lim = (g_avg), (5.0*m.log10(x_lim[1]/0.01) + g_avg)
    ax2 = ax1.twiny()
    ax2.set_xlim(new_lim)
    plt.xticks(x_tick, new_ticks)
    plt.xticks(fontsize=10)
    plt.xlabel(r'$g_{0}$', fontsize=10)
    plt.show()
    return 0

if __name__ == "__main__":
    plot_stuff(sys.argv[1:])