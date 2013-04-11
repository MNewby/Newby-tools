#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import gradient_descent as gd
import gauss_fit_enc as gfe

"""This program contains code for finding errors from R-squared functions.
Simply steps away from best fit parameters, finds R-squared, then fits a
normal distribution to the points.
Matthew Newby (RPI), Oct 16, 2010
"""

#Need to flip y values before evaluation!!!
def gaussian_function (x_in, params):
    exponent = -0.5*( (x_in - params[0]) / params[1])**2
    factor = 1.0 /(np.sqrt(2.0*m.pi*params[1]*params[1]))
    y_out = factor*np.exp(exponent)
    return y_out

#Need - function, params, step size
#Points is the number of points on each side of the fit parameter
#Steps should be an array the same size as best_params
def get_normal_errors(function, best_params, hist_in, sigma, steps, points=10):
    #replace hist_in with x,y
    print '#---Calculating errors'
    spread = 2*points + 1
    error_points = sc.zeros((len(best_params),spread),float)
    for i in range(len(best_params)):
        error_points[i,0] = best_params[i]
        for j in range(points):
            error_points[i,j+1] = best_params[i] - (j+1)*steps[i]
        for j in range(points):
            error_points[i,j+points+1] = best_params[i] + (j+1)*steps[i]
    #get R-squared values
    Rsquared_points = sc.zeros((error_points.shape), float)
    l, w = error_points.shape
    for i in range(l):
        for j in range(w):
            test_params = sc.zeros(len(best_params))
            for k in range(len(best_params)):  test_params[k] = best_params[k]
            test_params[i] = error_points[i,j]
            Rsquared_points[i,j] = R_squared_gauss(hist_in, test_params)
            """Rsquared_points[i,j] = R_squared(function, test_params,x,y,sigma)"""
    #Find gaussian fit to each parameter - maybe only fit sigma?
    error_out = []
    for i in range(len(best_params)):
        error_x = error_points[i]
        error_y = Rsquared_points[i]
        error_y = -1.0*(error_y / np.ma.max(error_y)) + 1.0
        error_sigma = sc.ones(len(error_x))
        error_params = sc.array([best_params[i],steps[i]])
        error_fit = gd.gradient_descent(gaussian_function, error_params,error_x,
                                        error_y, error_sigma)
        #error_fit = gfe.Gradient_Descent(sc.array([error_x, error_y]), error_params)
        print '#---Error for parameter #', i
        print '#-parameter:', error_x
        print '#-R-squared:', error_y
        print '#-Error fit mean:', error_fit[0], ', sigma:', error_fit[1]
        error_out.append(error_fit)
    return sc.array(error_out)

def R_squared(function, parameters, x, y, sigma):
    #print '-R squared params', parameters
    pieces = ( (y - function(x, parameters)) / sigma )**2
    return sc.sum(pieces)
    
#From gauss_fit_enc.py
def R_squared_gauss (hist_in, parameters):
    #This function returns the r-squared value of a data set, x,y, that is "points" long,
    #and the current gaussian parameter set, mu, sigma.  Note that this gaussian is double-sided.
    #amplitude is the height of the gaussian
    length, width = hist_in.shape
    mu, sigma_left, sigma_right, amplitude = parameters
    line_y = sc.zeros(length)
    #Create the gaussian from parameter sets
    for i in range(length):
        if hist_in[i,0] < mu:
            sigma = sigma_left
        else:
            sigma = sigma_right
        exponent = -0.5 * ( (hist_in[i,0] - mu) / sigma )**2
        stuff = m.exp(exponent)
        line_y[i] = stuff*amplitude
    #Find r-squared value
    r_squared = 0.0
    for i in range(length):
        if (line_y[i] < 1.0):
            denom = 1.0
        else:
            denom = line_y[i]
        diff = ((hist_in[i,1] - line_y[i])*(hist_in[i,1] - line_y[i]))/ (denom)
        r_squared = r_squared + diff
    return r_squared