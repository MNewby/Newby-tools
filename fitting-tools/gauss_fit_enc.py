#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import time

"""This program contains code for fitting gaussian distributions,
   and is nicely encapsulated for use.

Matthew Newby (RPI), June 3, 2010
"""

#Performs an MCMC fit to a double-sided gaussian data set.
#params_old is starting parameters.
def do_MCMC(hist_data, params_old, iterations=10000):
    MCMC_factor = 8.0
    randomSeed = int(time.time())
    np.random.seed(randomSeed)
    RR_old = R_squared_gauss (hist_data, params_old)
    RR_best = RR_old
    params_best = params_old
    print '-Starting MCMC for', iterations, 'iterations'
    for i in range(iterations):
        RR_new, params_new = MCMC_step(hist_data, params_old)
        if (RR_new < RR_best):
            RR_best = RR_new
            params_best = params_new
        if (RR_new < RR_old):
            RR_old = RR_new
            params_old = params_new
        else: 
            if ( (RR_old / (RR_new*MCMC_factor)) > np.random.uniform() ):
                RR_old = RR_new
                params_old = params_new
        if (i % 1000) == 0: 
            print i, params_old, RR_old        
    return RR_best, params_best

def Gradient_Descent(hist_in, st_parameters=[0.0,0.0,0.0,0.0], iterations=1000, precision=0.000000001):
    #Finds the best fit to a double-sided gaussian.
    #May produce bad results if Amp. is small; finds strange local min.
    #for manual starting point, beware of new params being better (use new_params instead):
    print '-Starting Gradient Descent...'
    print st_parameters
    old_params = [100.0, 100.0, 100.0, 100.0]
    new_params = st_parameters
    gradient = [1000.0, 1000.0, 1000.0, 1000.0]
    RR_old = R_squared_gauss(hist_in, old_params)
    RR_new = R_squared_gauss(hist_in, new_params)
    loop = 0
    overrun = 0
    scale = 0.01
    while (abs(RR_old-RR_new) > precision):  #MAY STILL MISS LAST STEP!???
        if (RR_new < RR_old):
            print '### I moved!  r-squared change:', abs(RR_old - RR_new)
            RR_old = RR_new
            for k in range(4):
                old_params[k] = new_params[k]
            scale = scale*1.5
        else:
            scale = scale*0.5
        loop = loop + 1
        if (loop > iterations):
            print '-Exited by passing iteration threshold'
            overrun = 1
            break
        gradient[0] = DR2_Dmu(hist_in, old_params)
        gradient[1], gradient[2] = DR2_Dsigma(hist_in, old_params)
        gradient[3] = DR2_Damp(hist_in, old_params)
        for i in range(4):
            new_params[i] = old_params[i] - scale*gradient[i]
        RR_new = R_squared_gauss(hist_in, new_params)
    if overrun == 0:
        print '-Gradient Descent successful, exited after ', loop, ' iterations'
    return old_params

#Takes in x, y data array and chooses gaussian parameters to start from.
def make_start_params(data):
    max_height = 0.0
    mu_n = 0
    length, width = data.shape
    for i in range(length):
        if (data[i,1] > max_height):
            max_height = data[i,1]
            mu_n = i
    mu = data[mu_n, 0]
    sigma_l, sigma_r = 1.0, 1.0
    amp = max_height
    start_params = [mu, sigma_l, sigma_r, amp]
    return start_params

def MCMC_step(hist_in, parameters, double_sided=1):
    #fits a gaussian curve to the data set, using MCMC
    new_params = [0.0, 0.0, 0.0, 0.0]
    x_range = (np.ma.max(hist_in[:,0]) - np.ma.min(hist_in[:,0]))
    y_range = (np.ma.max(hist_in[:,1]) - np.ma.min(hist_in[:,1]))
    new_params[0] = np.random.normal(parameters[0], (x_range/100.0))
    new_params[1] = np.random.normal(parameters[1], (x_range/100.0))
    if (double_sided == 1):
        new_params[2] = np.random.normal(parameters[2], (x_range/100.0))
    else:
        new_params[2] = new_params[1]
    new_params[3] = np.random.normal(parameters[3], (y_range/100.0))
    r_squared = R_squared_gauss(hist_in, new_params)
    return r_squared, new_params

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
        denom = poisson_errors(hist_in[i,1])
        #denom = poisson_errors(line_y[i])
        #if (line_y[i] < 1.0):
        #    denom = 1.0
        #else:
        #    denom = line_y[i]
        #denom = 1.0
        diff = ((hist_in[i,1] - line_y[i])*(hist_in[i,1] - line_y[i]))/ (denom)
        r_squared = r_squared + diff
    return r_squared

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

def DR2_Dmu(hist_in, params):
    #returns the gradient of R^2 with respect to mu at current mu
    #params[mu, sigma_l, sigma_r, amplitude] 
    mu_grad = 0.0000
    length, width = hist_in.shape
    for i in range(length):
        if (hist_in[i,0] < params[0]):
            sigma = params[1]
        else:
            sigma = params[2]
        part1 = (hist_in[i,0] - params[0])/(sigma*sigma)
        part2 = gaussian(hist_in[i,0], params[0], sigma, params[3])
        part3 = (hist_in[i,1] - part2)
        add_on = part1 * part2 * part3
        mu_grad = mu_grad + add_on
    return (mu_grad*(-2.0))


def DR2_Dsigma(hist_in, params):
    #returns the gradient of R^2 with respect to sigma at current sigma
    #params[mu, sigma_l, sigma_r, amplitude] 
    sigma_l_grad = 0.0000
    sigma_r_grad = 0.0000
    length, width = hist_in.shape
    for i in range(length):
        if (hist_in[i,0] < params[0]):
            sigma = params[1]
        else:
            sigma = params[2]
        part1 = ((hist_in[i,0] - params[0])**2)/(sigma*sigma*sigma)
        part2 = gaussian(hist_in[i,0], params[0], sigma, params[3])
        part3 = (hist_in[i,1] - part2)
        add_on = part1 * part2 * part3
        if (hist_in[i,0] < params[0]):
            sigma_l_grad = sigma_l_grad + add_on
        else:
            sigma_r_grad = sigma_r_grad + add_on
    return (sigma_l_grad*(-2.0)), (sigma_r_grad*(-2.0))
 
def DR2_Damp(hist_in, params):
    #returns the gradient of R^2 with respect to amplitude at current amp.
    #params[mu, sigma_l, sigma_r, amplitude] 
    amp_grad = 0.0000
    length, width = hist_in.shape
    for i in range(length):
        if (hist_in[i,0] < params[0]):
            sigma = params[1]
        else:
            sigma = params[2]
        part1 = m.exp(-0.5*( (hist_in[i,0] - params[0]) / sigma)**2)
        part2 = (hist_in[i,1] - gaussian(hist_in[i,0], params[0], sigma, params[3]) )
        add_on = part1 * part2
        amp_grad = amp_grad + add_on
    return (amp_grad*(-2.0))

def gaussian (x_in, mu, sigma, amplitude=1.0):
    exponent = -0.5*( (x_in - mu) / sigma)**2
    y_out = amplitude*m.exp(exponent)
    return y_out

def full_width_half_max (sigma_l, sigma_r=0):
    if (sigma_r == 0):
        sigma_r = sigma_l
    width_l = (2.35482/2.0)*m.fabs(sigma_l)
    width_r = (2.35482/2.0)*m.fabs(sigma_r)
    FWHM = width_l + width_r
    return FWHM

def poisson_errors(n, rigorous=False):
    #upper error bar
    lambda_up = n + np.sqrt(n+1.0) + 1.0
    #lower error bar
    #lambda_down = n*( (1.0 - (1.0/(9.0*n)) - (1.0/(3.0*np.sqrt(n))))**3 )
    return lambda_up

def make_data(real_mu=5.0, real_sigma_l=2.5, real_sigma_r=2.5, size = 30, error = 2.0, amplitude = 10.0):
    #Makes a gaussian with random errors for testing purposes
    randomSeed = int(time.time())
    np.random.seed(randomSeed)
    holder = np.random.rand(size)
    holder = ( (holder*error*2) - error)
    x = sc.zeros(size)
    y = sc.zeros(size)
    #Make x-range and y-range for data set
    length = real_sigma_r*4.0
    minimum = real_mu - length
    step = (2.0*length) / size
    for j in range(size):
        x[j] = minimum + (step*j)
        if (x[j] < real_mu):
            sigma = real_sigma_l
        else:
            sigma = real_sigma_r
        exponent = -0.5*( (x[j] - real_mu) / sigma)**2
        y[j] = ( m.exp(exponent)*amplitude) + holder[j]
    data_out = sc.zeros((size, 2))
    data_out[:,0] = x
    data_out[:,1] = y
    return data_out

#########################################################################################
#    if GD == 1:
#        params_best = Gradient_Descent(x, y, size, params_old)
#        RR_best = R_squared_gauss(x,y,size, params_best)
#########################################################################################







