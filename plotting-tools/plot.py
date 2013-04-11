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

args:  Filename, plot type:
    scatter, x column, y column
    hist, data column, bins

Matthew Newby (RPI), Dec 10, 2010
"""

filename = sys.argv[1]
plot_type = sys.argv[2]
suffix = filename[-4:]
if suffix == '.txt':
    data = f.read_data(filename)
elif suffix == '.csv':
    data = f.read_data(filename, ',') 

plt.figure(1)
'''Plot Type Selection'''
if plot_type == 'scatter':
    print "Plotting scatter plot"
    x_col = sys.argv[3]
    y_col = sys.argv[4]
    data_x = data[:,x_col]
    data_y = data[:,y_col]
    plt.scatter(data_x, data_y, 1, 'k', 'o')
elif plot_type == 'cmd':
    print "plotting simple cmd (color1 vs. color1-color2)"
    col_1 = sys.argv[3]
    col_2 = sys.argv[4]
    data_1 = data[:,col_1]
    data_2 = data[:,col_2]
    new_x = (data_1 - data_2)
    plt.scatter(new_x, data_1, 1, 'k', 'o')
    #plt.plot([-0.3, 0.6],[23.0, 23.0], 'r:')
    #limits = plt.ylim()
    #plt.ylim(limits[1], limits[0])
    plt.ylim(26.0, 15.0)
    plt.xlim(0.0, 1.0)
    plt.ylabel(r'$g_0$')
    plt.xlabel(r"$(g-r)_0$")
elif plot_type == 'hist':
    print "Plotting histogram"
    x_col = sys.argv[3]
    nbins = sys.argv[4]
    data_x = data[:,x_col]
    plt.hist(data_x, nbins)  #20
else:
    #if type(plot_type) != type(1):
    #    print "Second argument doesn't make sense, please double-check it"
    #    sys.exit(2)
    print 'Do not understand plot_type; assuming series plot'
    x_col = sys.argv[2]
    y_col = sys.argv[3]
    data_x = data[:,x_col]
    data_y = data[:,y_col]
    plt.plot(data_x, data_y, 'k:')

#plt.xlabel('Distance (kpc)')
#plt.ylabel('Saturation Magnitude, $g$')

"""Actual Plotting Code"""
#plt.errorbar(data_x, data_y, yerr=y_err, ecolor='k', fmt=None)
#plt.scatter(data_x, data_y, marker='+')
#plt.plot(line_x, line_y, 'r--')
#plt.ylim(-0.1, 0.5)
#plt.xlim(12.0, 24.0)
#plt.axis('scaled')  #Scales the plot to make it the axes proportional
#plt.plot(x,y)
plt.show()
plt.close('all')
print '#-Done'

def gaussian_function (x_in, params):
    exponent = -0.5*( (x_in - params[0]) / params[1])**2
    factor = 1.0 /(np.sqrt(2.0*m.pi*params[1]*params[1]))
    y_out = factor*np.exp(exponent)
    return y_out

def detection_sigmoid(x_in, params=[]):
    s = [0.9402, 1.6171, 23.5877]
    y_out = s[0] / (np.exp(s[1]*(x_in - s[2])) + 1.)
    return y_out

supported = "scatter, cmd, line, hess"
supported_inputs = " scatter or line - x:y\n cmd or hess - ymag:xcolor1:xcolor2\n"

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="infile",
                      help="file to be read in, default is first arg", default=sys.argv[1])
    parser.add_option("-t", "--type", dest="plot_type", default='scatter',
                      help="the type of plot to generate.  currently supported: {0}".format(supported))
    parser.add_option("-c", "--columns", dest="columns", default="0:1",
                      help="columns to plot. Inputs depend on plot type: {0}".format(supported_inputs))
    parser.add_option("-x", "--xlabel", dest="xlabel", default="x",
                      help="label for the x axis, accepts latex formatting?")
    parser.add_option("-y", "--ylabel", dest="ylabel", default="y",
                      help="label for the y axis, accepts latex formatting?")
    parser.add_option("--xlim", dest="xlim", default=None,
                      help="range for the x axis - low:high")
    parser.add_option("--ylim", dest="ylim", default=None,
                      help="range for the y axis - low:high")   
    options, extra_args = parser.parse_args(args)
    #MORE?
    #Process args