#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time

"""This program contains code for fitting gaussian distributions.
I should really make this more generic...

Matthew Newby (RPI), Jan 28, 2010
"""

randomSeed = int(time.time())
np.random.seed(randomSeed)
GD = 1
unknown_GD = 0  #not working yet
MCMC = 0
iterations = 10000
precision = .000000001  #gradient threshold for gradient descent to exit
MCMC_factor = 8.0   #Influences movement prob. for MCMC, higher->less likely
GDstart = [19.3911434611, 0.19440684351639326, 1.2056532633265824, 22.005153886721963]

def test():
    real_mu = 5.0
    real_sigma_l = 1.0
    real_sigma_r = 3.5
    size = 40
    error = 3.0
    amplitude = 40.0
    x_values, y_values = make_data(real_mu, real_sigma_l, real_sigma_r, size, error, amplitude)

    line_params = gauss_fit(x_values, y_values, size, 1,1,0,0.0,0.0)

    axis = (np.ma.max(x_values) - np.ma.min(x_values))/2.0
    x_line = sc.arange( (np.ma.min(x_values)-axis), (np.ma.max(x_values)+axis), 0.01)
    y_line = sc.zeros(len(x_line))
    for j in range(len(x_line)):
        if x_line[j] < line_params[0]:
            sigma = line_params[1]
        else:
            sigma = line_params[2]
        exponent = -0.5 * ( (x_line[j] - line_params[0]) / sigma )**2
        stuff = m.exp(exponent)
        y_line[j] = stuff*line_params[3]

    print 'best parameters:', line_params
    plt.figure(1)
    plt.scatter(x_values, y_values)
    plt.plot(x_line, y_line, 'g-')
    plt.show()
    return 'Done'


#put bool and floats in one nice package (list) each
def gauss_fit(x,y,size, fit_amp=0, double_sided=0, cut=0, low_cut=0.0, high_cut=100.0):
    #fits a gaussian to the data set.  fit_amp is 1 if you want to fit the amplitude,
    #double_sided is 1 if your gaussian is double sided.  Cut is 1 if you want to fit to
    #a limited x range; low_cut and high_cut are the range in which to fit.  CUT NOT WORKING YET
    
    #Create starting values:
    maximum = 0.0
    mu_n = 0
    for j in range(size):
       if (y[j] > maximum):
           maximum = y[j]
           mu_n = j
    mu = x[mu_n]
    sigma_l = 1.0
    sigma_r = 1.0
    amplitude = maximum
    params_old = [mu, sigma_l, sigma_r, amplitude]

    #time for a least-squares fit:
    RR_old = R_squared_gauss(x,y,size, params_old)
    RR_best = RR_old
    params_best = params_old

    #Monte-Carlo Markov Chain - move into a MCMC function
    if MCMC == 1:
        print '-Starting MCMC for', iterations, 'iterations'
        print '-With seed#:', randomSeed, 'and MC factor', MCMC_factor, '...'
        for k in range(iterations):
           RR_new, params_new = MCMC_step(x,y,size, params_old, double_sided)
           if (RR_new < RR_best):
               RR_best = RR_new
               params_best = params_new
           if (RR_new < RR_old):
               RR_old = RR_new
               params_old = params_new
               #print k, params_old, RR_old
           else: 
               if ( (RR_old / (RR_new*MCMC_factor)) > np.random.uniform() ):
                   RR_old = RR_new
                   params_old = params_new
                   #print k, params_old, RR_old
           if (k % 100) == 0: 
               print k, params_old, RR_old

    if GD == 1:
        params_best = Gradient_Descent(x, y, size, params_old)
        RR_best = R_squared_gauss(x,y,size, params_best)

    #Parameters are not iterating - gradient is too big!!!
    if unknown_GD == 1:
        loop = 0
        steps = [0.0, 0.0, 0.0, 0.0]
        steps[0] = (np.ma.max(x) - np.ma.min(x)) / 100.0
        steps[1] = (np.ma.max(x) - np.ma.min(x)) / 100.0
        steps[2] = (np.ma.max(x) - np.ma.min(x)) / 100.0
        steps[3] = (np.ma.max(y) - np.ma.min(y)) / 100.0
        print 'steps: ', steps
        scale = [0.001, 0.001, 0.001, 0.001]
        while (np.ma.max(scale) > scale_threshold):  
            loop = loop + 1
            if (loop > iterations):
                print '-Exited by passing iteration threshold'
                break
            RR_new, params_new = unknown_GD_step(x,y,size, params_old, steps, scale)
            print 'RR: ', RR_old, RR_new
            if (RR_new < RR_old):
                RR_old = RR_new
                params_old = params_new
            #print loop, params_old, RR_old, scale

        RR_best = RR_old
        params_best = params_old
        
    print '-best r-squared:', RR_best
    return params_best

def Gradient_Descent(x, y, size, old_params):
    #Finds the best fit to a double-sided gaussian.
    #RR_new is initially set to new_params[0s] so that difference will be large.
    #May produce bad results if Amp. is small; finds strange local min.
    #for manual starting point, beware of new params being better (use new_params instead):
    #old_params = [50.00, 1.0, 1.0, 80.00]
    print '-Starting Gradient Descent...'
    print old_params
    gradient = [1000.0, 1000.0, 1000.0, 1000.0]
    #use the first line below for an intentional "bad" start (default); 2nd line to choose.
    #new_params = [1000.0, 1000.0, 1000.0, 1000.0]
    new_params = GDstart
    RR_old = R_squared_gauss(x, y, size, old_params)
    RR_new = R_squared_gauss(x, y, size, new_params)
    loop = 0
    overrun = 0
    scale = 0.01
    while (abs(RR_old-RR_new) > precision):  #MAY STILL MISS LAST STEP!???
    #while (scale > 0.0000000000000001):
        #print 'old params, r-squared:', old_params, RR_old
        #print '-gradient, scale: ', gradient, scale
        #print '-new params, r-squared: ', new_params, RR_new
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
        gradient[0] = DR2_Dmu(x, y, old_params)
        gradient[1], gradient[2] = DR2_Dsigma(x, y, old_params)
        gradient[3] = DR2_Damp(x, y, old_params)
        for k in range(4):
            new_params[k] = old_params[k] - scale*gradient[k]
        RR_new = R_squared_gauss(x, y, size, new_params)
        #if (loop % 10) == 0: 
        #    print loop, old_params, gradient, scale, RR_old
    if overrun == 0:
        print '-Gradient Descent successful, exited after ', loop, ' iterations'
    return old_params


def unknown_GD_step(x, y, size, old_params, steps, scale):
    #CRAP - needs cleaning to be useful; use only for unknown gradients.
    #Find the next gradient descent position, inputs are starting parameters, 
    #parameter step sizes, and the gradient scaling factor.
    new_params = [0.0, 0.0, 0.0, 0.0]
    test_params_plus = [0.0, 0.0, 0.0, 0.0]
    test_params_minus = [0.0, 0.0, 0.0, 0.0]
    old_RR = R_squared_gauss(x, y, size, old_params)
    for j in range(len(old_params)):
        for k in range(len(old_params)):
            test_params_plus[k] = old_params[k]
            test_params_minus[k] = old_params[k]
        #print 'p-m: ', test_params_plus[j], test_params_minus[j]
        plus_value = old_params[j] + steps[j]
        minus_value = old_params[j] - steps[j]
        #print 'values ', plus_value, minus_value
        test_params_plus[j] = plus_value
        test_params_minus[j] = minus_value
        #print 'p-m: ', test_params_plus[j], test_params_minus[j], steps[j]
        plus = R_squared_gauss(x, y, size, test_params_plus)
        minus = R_squared_gauss(x, y, size, test_params_minus)
        gradient = (plus - minus)/(2.0*steps[j])
        print 'p,m, grad: ', plus, minus, gradient
        new_params[j] = old_params[j] - (scale[j]*gradient)
        better_test = R_squared_gauss(x, y, size, new_params)
        #print gradient, scale[j], better_test
        if (better_test < old_RR):
            scale[j] = 1.2*scale[j]
        else:
            scale[j] = 0.8*scale[j]
            new_params[j] = old_params[j]

    r_squared = R_squared_gauss(x, y, size, new_params)
    return r_squared, new_params


def MCMC_step(x_data, y_data, size, parameters, double_sided):
    #fits a gaussian curve to the data set, using MCMC
    st_dev = [0.0, 0.0, 0.0, 0.0]
    new_params = [0.0, 0.0, 0.0, 0.0]
    st_dev[0] = (np.ma.max(x_data) - np.ma.min(x_data)) / 100.0
    st_dev[1] = (np.ma.max(x_data) - np.ma.min(x_data)) / 100.0
    st_dev[2] = (np.ma.max(x_data) - np.ma.min(x_data)) / 100.0
    st_dev[3] = (np.ma.max(y_data) - np.ma.min(y_data)) / 100.0
    new_params[0] = np.random.normal(parameters[0], st_dev[0])
    new_params[1] = np.random.normal(parameters[1], st_dev[1])
    if (double_sided == 1):
        new_params[2] = np.random.normal(parameters[2], st_dev[2])
    else:
        new_params[2] = new_params[1]
    new_params[3] = np.random.normal(parameters[3], st_dev[3])
    r_squared = R_squared_gauss(x_data, y_data, size, new_params)
    return r_squared, new_params


def R_squared_gauss (data_x, data_y, points, parameters):
    #This function returns the r-squared value of a data set, x,y, that is "points" long,
    #and the current gaussian parameter set, mu, sigma.  Note that this gaussian is double-sided.
    #amplitude is the height or the gaussian

    mu, sigma_left, sigma_right, amplitude = parameters
    line_y = sc.zeros(points)
    
    #Create the gaussian from parameter sets
    for j in range(points):
        if data_x[j] < mu:
            sigma = sigma_left
        else:
            sigma = sigma_right
        exponent = -0.5 * ( (data_x[j] - mu) / sigma )**2
        stuff = m.exp(exponent)
        line_y[j] = stuff*amplitude

    #Find r-squared value
    r_squared = 0.0
    for k in range(points):
        if (line_y[k] < 1.0):
            denom = 1.0
        else:
            denom = line_y[k]
        diff = ((data_y[k] - line_y[k])*(data_y[k] - line_y[k]))/ (denom)
        r_squared = r_squared + diff

    return r_squared


def make_data(real_mu=5.0, real_sigma_l=2.5, real_sigma_r=2.5, size = 30, error = 1.0, amplitude = 4.0):
    #Makes a gaussian with random errors for testing purposes
    x = sc.zeros(size)
    y = sc.zeros(size)
    holder = np.random.rand(size)
    holder = ( (holder*error*2) - error)

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

    #plt.figure(1)
    #plt.scatter(x,y)
    #plt.show()
    return x, y


def full_width_half_max (sigma_l, sigma_r=0):
    if (sigma_r == 0):
        sigma_r = sigma_l
    width_l = (2.35482/2.0)*m.fabs(sigma_l)
    width_r = (2.35482/2.0)*m.fabs(sigma_r)
    FWHM = width_l + width_r
    return FWHM


def gaussian (x_in, mu, sigma, amplitude=1.0):
    exponent = -0.5*( (x_in - mu) / sigma)**2
    y_out = amplitude*m.exp(exponent)
    return y_out


def DR2_Dmu(x, y, params):
    #returns the gradient of R^2 with respect to mu at current mu
    #params[mu, sigma_l, sigma_r, amplitude] 
    mu_grad = 0.0000
    for j in range(len(x)):
        if (x[j] < params[0]):
            sigma = params[1]
        else:
            sigma = params[2]
        part1 = (x[j] - params[0])/(sigma*sigma)
        part2 = gaussian(x[j], params[0], sigma, params[3])
        part3 = (y[j] - part2)
        add_on = part1 * part2 * part3
        mu_grad = mu_grad + add_on
    return (mu_grad*(-2.0))


def DR2_Dsigma(x, y, params):
    #returns the gradient of R^2 with respect to sigma at current sigma
    #params[mu, sigma_l, sigma_r, amplitude] 
    sigma_l_grad = 0.0000
    sigma_r_grad = 0.0000
    for j in range(len(x)):
        if (x[j] < params[0]):
            sigma = params[1]
        else:
            sigma = params[2]
        part1 = ((x[j] - params[0])**2)/(sigma*sigma*sigma)
        part2 = gaussian(x[j], params[0], sigma, params[3])
        part3 = (y[j] - part2)
        add_on = part1 * part2 * part3
        if (x[j] < params[0]):
            sigma_l_grad = sigma_l_grad + add_on
        else:
            sigma_r_grad = sigma_r_grad + add_on
    return (sigma_l_grad*(-2.0)), (sigma_r_grad*(-2.0))
 

def DR2_Damp(x, y, params):
    #returns the gradient of R^2 with respect to amplitude at current amp.
    #params[mu, sigma_l, sigma_r, amplitude] 
    amp_grad = 0.0000
    for j in range(len(x)):
        if (x[j] < params[0]):
            sigma = params[1]
        else:
            sigma = params[2]
        part1 = m.exp(-0.5*( (x[j] - params[0]) / sigma)**2)
        part2 = (y[j] - gaussian(x[j], params[0], sigma, params[3]) )
        add_on = part1 * part2
        amp_grad = amp_grad + add_on
    return (amp_grad*(-2.0))


