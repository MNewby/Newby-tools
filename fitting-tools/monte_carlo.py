#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import time

"""This program contains code for fitting a generic function to data.
Needs the function and the data as imputs.

Matthew Newby (RPI), Sept 3, 2010
"""

def template_function(x, params):
    return (x*params[1] + params[0])
    
def do_MCMC(function, init_params, step_sizes, x, y, sigma, name, number_steps=10000,
            save=1):
    #name is a list of strings of this form: run name, param name 1, param name 2, ...
    np.random.seed( int(time.time()) )
    moves = 0
    positions = []
    best_RR = 10000000.0
    best_params = sc.zeros(len(init_params))
    current_params = sc.zeros(len(init_params))
    for i in range(len(init_params)):  current_params[i] = init_params[i]
    new_params = sc.zeros(len(init_params), float)
    for i in range(number_steps):
        #Record position
        positions.append(current_params.tolist())
        #Take a step
        for j in range(len(init_params)):
            new_params[j] = np.random.normal(current_params[j], step_sizes[j])
        #Decide whether to move or not
        current_RR = R_squared(function, current_params, x, y, sigma)
        new_RR = R_squared(function, new_params, x, y, sigma)
        compare = (current_RR / (new_RR*4.0))
        #if (np.random.uniform() < (np.exp(-1.0*new_RR) / np.exp(-1.0*current_RR) ) ):
        if ( (new_RR < current_RR) or (np.random.uniform() < compare) ):
            moves = moves + 1
            for j in range(len(init_params)):
                current_params[j] = new_params[j]
        """#Just move, dammit
        new_RR = R_squared(function, new_params, x, y, sigma)
        for j in range(len(init_params)):
                current_params[j] = new_params[j]
        moves = moves + 1
        #Remove from 'dammit' to here to restore old functionality"""
        #record best fitness
        if (new_RR < best_RR):
            best_RR = new_RR
            for j in range(len(init_params)):
                best_params[j] = new_params[j]
        #if (i % 1000) == 0:
        #    print 'At iteration', i, current_params
    #make histogram
    centers, deviations = plot_MCMC_hist(positions, name, save)
    print '#---Mean parameters:', centers, 'Parameter deviations:', deviations
    print '#---Best parameters:', best_params, best_RR
    print '#---After:', moves, ' moves out of ', number_steps,'iterations'
    #plot fit with function --- with both best and means?
    #if (plot_function(x, y, function, params, x_err=None, y_err=sigma,
    #              title=name[0], save=0, save_file=(name[0]+'_func_plot.ps') ) ==1):
    #    print '#---data and MCMC function fit successfully plotted'
    return centers, deviations, best_params

def R_squared(function, parameters, x, y, sigma):
    #print '-R squared params', parameters
    pieces = ( (y - function(x, parameters)) / sigma )**2
    return sc.sum(pieces)

def plot_MCMC_hist(positions, name, save=1, warmup=4000):
    #name is a list of strings of this form: run name, param name 1, param name 2, ...
    data = np.array(positions[warmup:], dtype=float)
    means, st_devs = sc.mean(data, 0), sc.std(data, 0)
    nbins = 20 #(len(positions) / 10)
    subplot_rows = int(len(positions[0])+0.5)
    plt.figure(1)
    plt.title(('Parameter plots of '+name[0]))
    #Create a subplot for each parameter
    for i in range(len(positions[0])):
        x_step = ( np.ma.max(data[:,i])-np.ma.min(data[:,i]) ) / 100.0
        line_x = sc.arange(np.ma.min(data[:,i]), np.ma.max(data[:,i]), x_step)
        line_y = ( ( 1/(np.sqrt(2*m.pi)*st_devs[i]) )*
            (np.exp(-0.5*( (line_x - means[i])/st_devs[i])**2) ) )
        print '#---Making subplot:', subplot_rows, 2, i+1
        plt.subplot(subplot_rows, 2, (i+1))
        plt.hist(data[:,i], nbins, normed=True)
        plt.plot(line_x, line_y)
        plt.xlabel(name[i+1])
    if save == 1:
        save_file = name[0] + '_MCMChist_.ps'    
        plt.savefig(save_file, papertype='letter')
        print '#---MCMC hist plot saved as', save_file
    elif save == 0:
        pass
    else:
        plt.show()
    plt.close('all')
    print '#---Succesfully Plotted ', name[0], 'param plots'
    return means, st_devs


"""------------WORKING LINE----------------

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
        if (line_y[i] < 1.0):
            denom = 1.0
        else:
            denom = line_y[i]
        diff = ((hist_in[i,1] - line_y[i])*(hist_in[i,1] - line_y[i]))/ (denom)
        r_squared = r_squared + diff
    return r_squared

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
"""







