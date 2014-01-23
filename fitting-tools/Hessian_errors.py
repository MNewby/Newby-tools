#! /usr/bin/env python

import scipy as sc
import numpy as np

"""This program contains code for finding errors from R-squared functions.
Matthew Newby (RPI), Oct 16, 2010
http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node8.html
N_data = number of data points in original data
"""

def get_hessian_errors(function, best_params, x, y, sigma, step, like=0):
    length = len(best_params)
    hessian = sc.zeros((length, length),float)
    base_params = sc.zeros(length, float)
    for i in range(length): base_params[i] = best_params[i]
    test_params = sc.zeros(length, float)
    #Build Hessian
    for i in range(length):
        for j in range(length):
            #H_1[i,j] = L(Q[j]+h[j], Q[i]+h[i])
            for k in range(length): test_params[k]=base_params[k]
            test_params[j] = test_params[j] + step[j]
            test_params[i] = test_params[i] + step[i]
            H_1 = R_squared(function, test_params, x, y, sigma)
            if like==1:  H_1 = np.exp(-H_1/2.0)
            #H_2[i,j] = L(Q[j]-h[j], Q[i]+h[i])
            for k in range(length): test_params[k]=base_params[k]
            test_params[j] = test_params[j] - step[j]
            test_params[i] = test_params[i] + step[i]
            H_2 = R_squared(function, test_params, x, y, sigma)
            if like==1:  H_2 = np.exp(-H_2/2.0)
            #H_3[i,j] = L(Q[j]+h[j], Q[i]-h[i])
            for k in range(length): test_params[k]=base_params[k]
            test_params[j] = test_params[j] + step[j]
            test_params[i] = test_params[i] - step[i]
            H_3 = R_squared(function, test_params, x, y, sigma)
            if like==1:  H_3 = np.exp(-H_3/2.0)
            #H_4[i,j] = L(Q[j]-h[j], Q[i]-h[i])
            for k in range(length): test_params[k]=base_params[k]
            test_params[j] = test_params[j] - step[j]
            test_params[i] = test_params[i] - step[i]
            H_4 = R_squared(function, test_params, x, y, sigma)
            if like==1:  H_4 = np.exp(-H_4/2.0)
            #H[i,j] = (H_1[i,j] - H_2[i,j]-H_3[i,j]+H_4[i,j]) / (4*h[i]*h[j])
            #h[k] is the step size for likelihood determination; use same as gradient?
            hessian[i,j] = (H_1 - H_2 - H_3 + H_4) / (4.0*step[i]*step[j])
    #check that it is symmetric
    for i in range(length):
        for j in range(length):
            if (hessian[i,j] != hessian[j,i]): print '!!! Hessian not symmetric!'
    #Convert to matrix type and invert
    hessian_matrix = np.matrix(hessian)
    print '#--Hessian Matrix:'
    print hessian_matrix
    hessian_inverse = hessian_matrix.I
    #read off diagonals and calculate errors
    errors = sc.zeros(length, float)
    for i in range (length):
        if (hessian_inverse[i,i] < 0.0):  print '!!!Hessian error', i, 'below zero!!!'
        errors[i] = np.sqrt(2.0*abs(hessian_inverse[i,i]))  #*(1.0/len(x))
        #errors[i] = np.sqrt(2.0*abs(hessian_inverse[i,i]))  #*(1.0/len(x))
    print '#---Hessian Errors:', errors
    return errors

def R_squared(function, parameters, x, y, sigma):
    pieces = ( (y - function(x, parameters)) / sigma )**2
    return sc.sum(pieces)
    
def reduced_chi_squared(chi_squared, N_data, n_params):
    return ( chi_squared/(N_data-n_params) )
