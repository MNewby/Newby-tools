#! /usr/bin/env python  #'Bang' line - modify as needed

import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import sys

"""This is a sample script for plotting data using matplotlib, which requires
scipy and numpy to run.

args:  Filename, x column, y column

Matthew Newby (RPI), Dec 10, 2010
"""

formats = ['b:', 'g-', 'r:', 'k-', 'c-', 'm-', 'b-', 'r:', 'g-', 'k-']

def get_2gauss_y (x, parameters):
    #Creates Y values for given x values for a double-sided gaussian function.
    length = len(x)
    mu, sigma_left, sigma_right, amplitude = parameters
    line_y = sc.zeros(length)
    #Create the gaussian from parameter sets
    for i in range(length):
        if x[i] < mu:  sigma = sigma_left
        else:  sigma = sigma_right
        exponent = -0.5 * ( (x[i] - mu) / sigma )**2
        stuff = m.exp(exponent)
        line_y[i] = stuff*amplitude
    return line_y

def plot_gauss_series(parameters):
    plt.figure(1)
    x = sc.arange(0.0, 10.0, 0.1)
    for i in range(len(parameters[:,0])):
        params = [parameters[i,2], parameters[i,3], parameters[i,4], parameters[i,5]]
        y = get_2gauss_y(x, params)
        plt.errorbar(x, y, yerr=None, fmt='k-')#, ecolor='g', marker='D', markersize=3)
    plt.xlabel(r'$M_{g}$')
    plt.ylabel('counts')
    plt.xlim(0.0, 10.0)
    plt.show()
    plt.close('all')
    print '#-Done'    
    return 0

if __name__ == '__main__':
    # "NGC_6205_series.txt"
    filename = sys.argv[1]
    start_line = (int(sys.argv[2]) - 1)
    end_line = int(sys.argv[3])
    data = f.read_data(filename)
    plot_gauss_series(data[start_line:end_line,:])





