#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt

"""This program plots data and a function on the same plot

Matthew Newby (RPI), Sept 8, 2010
"""

def plot_function(data_x, data_y, function, params, x_err=None, y_err=None,
                  title=['', 'x', 'y'], save=1, save_name='func_plot'):
    #Get ranges and expand them a bit
    x_min, x_max = np.ma.min(data_x), np.ma.max(data_x)
    run = (x_max - x_min)
    #x_min, x_max = (x_min - run*0.1), (x_max + run*0.1)
    #print x_min, x_max
    #step = (x_max - x_min)/1000.0
    #print step
    line_x = np.arange(x_min, x_max, 0.01)
    line_y = function(line_x, params)
    #print line_x, '\n', line_y
    plt.figure(1)
    plt.errorbar(data_x, data_y, yerr=y_err, xerr=x_err, ecolor='b', fmt=None)
    plt.scatter(data_x, data_y, s=1, marker='o', c='black')
    plt.plot((line_x), (line_y), c='g')
    #plt.title(title[0])
    plt.xlabel(title[1], fontsize=10)
    plt.ylabel(title[2], fontsize=10)
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
    import functions as func
    import files as f
    data = f.read_data('test_data_out.txt')
    parameters = [ 5.0, 0.0, 2.0]
    plot_function(data[:,0], data[:,1], func.Nth_order_polynomial, parameters, y_err=data[:,2],
                  save=0)