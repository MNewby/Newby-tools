#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import files as f
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt


"""This program plots data and a function on the same plot

Matthew Newby (RPI), Sept 8, 2010
"""

def sigmoid_function(x, params):
        #Looks good, takes three parameters
        data = ( params[0]/(1 + np.exp(-1.0*(x-params[1])) ) ) + params[2]
        return data

def plot_function(data_x, data_y, function, params, x_err=None, y_err=None,
                  title=['', 'x', 'y'], save=1, save_name='func_plot'):
    #Get ranges and expand them a bit
    x_min, x_max = np.ma.min(data_x), np.ma.max(data_x)
    run = (x_max - x_min)
    x_min, x_max = (x_min - run*0.25), (x_max + run*0.25)
    #print x_min, x_max
    #step = (x_max - x_min)/1000.0
    #print step
    line_x = np.arange(x_min, x_max, 0.1)
    line_y = function(line_x, params)
    #print line_x, '\n', line_y
    plt.figure(1)
    plt.errorbar(data_x, data_y, yerr=y_err, xerr=x_err, ecolor='b', fmt=None)
    plt.scatter(data_x, data_y, c='blue')
    plt.plot((line_x), (line_y), c='g')
    plt.title(title[0])
    plt.xlabel(title[1])
    plt.ylabel(title[2])
    if save == 1:
        save_file = save_name+'_pltfunc.ps'
        plt.savefig(save_file, papertype='letter')
        print '#---function and fit plot saved as', save_file
    else:
        plt.show()
    plt.close('all')
    print '#---Done with function and data plot'
    return 1

if __name__ == "__main__":
    data = f.read_data('sigr_fit_best_results.txt')
    parameters = [2.61973731,  2.66051657, -1.48711101]
    plot_function(
        data[:,0], data[:,4], sigmoid_function, parameters, y_err=data[:,8],
                  save_name = 'new')