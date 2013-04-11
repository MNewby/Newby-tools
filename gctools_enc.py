#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import gauss_fit_enc as ga
import time

"""This program contains code for manipulating .csv files, 
esp. for Globular Cluster analysis.
  -make_cut:  cuts the data set with a circular cut.
  -binned_background:  bins a 1D data set, and subracts background, if given
  -convolve: convolves a data set
  -select_stars: keeps only star in a certain range.

Matthew Newby (RPI), Jan 14, 2009
"""

#Definitions:::

def make_cut(data_in, center_ra, center_dec, inner_cut, outer_cut):
    #Takes in a 4-column data set (ra, dec, ...) and keeps only star within 
    #a circular cut. Cuts a ring out of the data, cetered on the gc center 
    #find the data points that belong in the new set    
    length, width = data_in.shape
    holder = sc.zeros(length, int)
    in_set = 0
    for i in range(length):
        d = m.sqrt( (data_in[i,0] - center_ra)**2 + (data_in[i,1] - center_dec)**2 ) 
        if (d > inner_cut) and (d < outer_cut):
            holder[i] = 1
            in_set = in_set + 1
    outliers = (length - in_set)
    print '#-Removed', outliers, 'from data set'
    #Make the new data set
    data_new = sc.zeros((in_set, width))
    j = 0
    for i in range(length):
        if holder[i] == 1:
            data_new[j,:] = data_in[i,:]
            j = j + 1
    if j != in_set:
        print '!!! ERROR: not all cut data written properly!'
    return data_new

def sigmoid_error(x, modulus=None):
    s = [0.9402, 1.6171, 23.5877]
    if modulus != None:
        s[2] = s[2] - modulus
    detection_efficiency = s[0] / (np.exp(s[1]*(x - s[2])) + 1.)
    return detection_efficiency

def parabola_error(x, modulus=None):
    params = [ -0.14267811, 20.92411611]
    if modulus != None:
        params[1] = params[1] - modulus
    if x < params[1]:
        detection_efficiency = 1.0
    elif (x > params[1] + np.sqrt(-1.0/params[0]) ):
        detection_efficiency = 0.0
    else:
        detection_efficiency = 1.0 + params[0]*((x - params[1])**2)
    return detection_efficiency

def binned_background(data_x, back_x=[], bin_size=0.1, in_area=1.0, back_area=1.0,
    cut_hist=0, low_cut=0.0, high_cut=10.0, threshold=None, modulus=None):
    #bins the data and the background x-data, then subtracts 
    #the background from the data (if a background exists)
    #Makes the range go from the minimum of both data sets to the maximum
    if back_x == []:
        x_min = np.ma.min(data_x)
        x_max = np.ma.max(data_x)
    else:
        if (np.ma.min(data_x) <= np.ma.min(back_x)): x_min = np.ma.min(data_x)
        else: x_min = np.ma.min(back_x)
        if (np.ma.max(data_x) >= np.ma.max(back_x)): x_max = np.ma.max(data_x)
        else: x_max = np.ma.max(back_x)
    nbins = int( ((x_max + bin_size) - x_min)/bin_size )  #add a bin to compensate for rounding
    print '#-', nbins, 'bins of size', bin_size, 'magnitudes, and a threshold value of', threshold
    #Creates the bins
    bin_centers = sc.zeros(nbins)
    bin_counts = sc.zeros(nbins)
    for i in range(nbins):
        bin_centers[i] = x_min + (i*bin_size) + (bin_size*0.5)
    #Bin data, and subtract background; no bin counts below zero
    u = 0  #Tracks the number of binned points from dataset
    v = 0  #Tracks the number of "unbinned" points, taken from background
    w = 0.0  #Tracks number of points removed due to hist_cut selection
    #bins all values between bin max and bin min, including points at bin min
    for i in range(nbins):
        bin_min = x_min + (i*bin_size)
        bin_max = x_min + ((i+1)*bin_size)
        for j in range(len(data_x)):
            if (bin_min <= data_x[j] < bin_max):
                bin_counts[i] = bin_counts[i] + 1.0
                u = u + 1  #bug test
        #Compensating for cluster rolloff
        if modulus != None:
            scale_denom = parabola_error(bin_centers[i], modulus)
            if ( round(scale_denom, 2) == 0.00):
                bin_counts[i] = 0.0
            else:
                bin_counts[i] = bin_counts[i] / scale_denom
    #Bins the last data point this inadvertantly adds a data point; not needed.
    '''for i in range(len(data_x)):
        if (data_x[i] == x_max):  
            bin_counts[-1] = bin_counts[-1] + 1.0
            u = u + 1  #tracking'''    
    #Removes background from bins
    if back_x != []:
        for i in range(nbins):
            bin_min = x_min + (i*bin_size)
            bin_max = x_min + ((i+1)*bin_size)
            back_bin = 0.0
            for j in range(len(back_x)):
                if (bin_min <= back_x[j] < bin_max):
                    back_bin = back_bin + 1.0
                    v = v + 1  #tracking
            #Compensating for background rolloff
            if modulus != None:
                scale_denom = sigmoid_error(bin_centers[i], modulus)
                if ( round(scale_denom, 2) == 0.00):
                    back_bin = 0.0
                else:
                    back_bin = back_bin / scale_denom
            #Compensating for area differences
            back_bin = back_bin*(in_area/back_area)
            bin_counts[i] = bin_counts[i] - back_bin
            if (bin_counts[i] < 0):
                bin_counts[i] = 0.00
        #Removes the last data point, if necessary 
        for i in range(len(back_x)):
            if (back_x[i] == x_max):  
                bin_counts[-1] = bin_counts[-1] - (in_area/back_area)
                v = v + 1  #tracking
            if (bin_counts[-1] < 0):
                bin_counts[-1] = 0.00
    #Cuts the data set, if selected
    if (cut_hist == 1):
        print '#-Using only data between M=', low_cut, 'and', high_cut
        for i in range(nbins):
            if (bin_centers[i] < low_cut):
                w = w + bin_counts[i]  #tracking
                bin_counts[i] = 0.0
            if (bin_centers[i] > high_cut):
                w = w + bin_counts[i]  #tracking
                bin_counts[i] = 0.0
    #Removes all bins below the threshold value, if a value was set
    if (threshold != None):
        below_thresh = 0
        for i in range(nbins):
            if (bin_counts[i] < threshold):
                bin_counts[i] = 0.0
                below_thresh = below_thresh + 1
    print '#-Total Data in:', len(data_x), '; total points binned:', u
    if (back_x != []):
        print '#-Total Background in:', len(back_x), '; total background removed:', v
    if (cut_hist == 1):
        print '#-Total number of points removed due to hist cuts:', w
    if (threshold != None):
        print '#-Total number of bins below threshold value:', below_thresh
    data_out = sc.zeros((nbins, 2))
    data_out[:,0] = bin_centers
    data_out[:,1] = bin_counts
    return data_out


#Takes a 4 column data set in, with g in the 3rd column and r in the 4th column.
#Convolves r and g to the SDSS errors at a new distance, 'con_dist'
#Returns a 4 column data set, with columns 1 & 2 unchanged, and g&r convolved.
#Possibly make this more generic, so g and r can be any column? 
def convolve(data_in, real_dist, con_dist, randomSeed=int(time.time())):
    np.random.seed(randomSeed)
    a_g = 0.0
    b_g = 0.790391
    c_g = -19.8928
    a_r = 0.0
    b_r = 0.766309
    c_r = -19.0334
    a_u = 0.02076735
    b_u = 0.81309147
    c_u = -19.61533299
    diff_mag = 5.0*(np.log10(con_dist / real_dist))
    print '#-Convolving data with random seed =', randomSeed
    print '#---from ', real_dist, 'kpc to', con_dist, 'kpc'
    length, width = data_in.shape
    g = sc.zeros(length)
    r = sc.zeros(length)
    u = sc.zeros(length)
    for i in range(length):
        sigma_g = np.sqrt((np.exp(b_g*(data_in[i,2]+diff_mag))**2 - np.exp(b_g*data_in[i,2])**2))*np.exp(c_g)
        sigma_r = np.sqrt((np.exp(b_r*(data_in[i,3]+diff_mag))**2 - np.exp(b_r*data_in[i,3])**2))*np.exp(c_r)
        sigma_u = np.sqrt(( (a_u + np.exp(b_u*(data_in[i,4]+diff_mag)+c_u))**2 - (a_u + np.exp((b_u*data_in[i,4]) + c_u))**2))
        if (sigma_g > 0.0):
            g[i] = np.random.normal(data_in[i,2], sigma_g)
        else:  
            g[i] = data_in[i,2]
        if (sigma_r > 0.0):
            r[i] = np.random.normal(data_in[i,3], sigma_r)
        else:  
            r[i] = data_in[i,3]
        if (sigma_u > 0.0):
            u[i] = np.random.normal(data_in[i,4], sigma_u)
        else:  
            u[i] = data_in[i,4]
    m = sc.zeros((length,5))
    for i in range(length):
        m[i,0] = data_in[i,0]
        m[i,1] = data_in[i,1]
        m[i,2] = g[i]
        m[i,3] = r[i]
        m[i,4] = u[i]
    return m

#Takes in a 4-column data file, with g and r as the last two columns,
# and returns only the stars in the given g-r range.
#could generalize to generic colors/column numbers.
def select_stars(data_in, low_limit = 0.0, high_limit = 0.3, ug_limit = 0.4):
    length, width = data_in.shape
    data_index = sc.zeros(length, int)    
    for i in range(length):
        if (data_in[i,2] - data_in[i,3]) >= low_limit:
            if (data_in[i,2] - data_in[i,3]) <= high_limit:
                data_index[i] = 1
        if (data_in[i,4] - data_in[i,2]) < ug_limit:
            data_index[i] = 0
    s = sc.sum(data_index)
    print '#-Stars in data set:', length, 'Stars in f-turnoff:', s
    new_data = sc.zeros((s,width))
    j = 0
    for i in range(length):
        if data_index[i] == 1:
            new_data[j,:] = data_in[i,:]
            j = j + 1
    if (j != s):  print '!!! stars out not equal to f-turnoff:', j, s
    return new_data

def select_stars_with_parabola(data_in, new_dist, randomSeed=int(time.time())):
    np.random.seed(randomSeed)
    keepers, sig_boot = 0, 0
    
    return 0

def select_stars_with_sigmoid(data_in, new_dist,
                              low_limit = 0.0, high_limit = 0.3, ug_limit = 0.4):
    np.random.seed(int(time.time()))
    s = [0.9402, 1.6171, 23.5877]
    length, width = data_in.shape
    data_index = sc.zeros(length, int)    
    keepers = 0
    sig_boot = 0
    for i in range(length):
        if (data_in[i,2] - data_in[i,3]) >= low_limit:
            if (data_in[i,2] - data_in[i,3]) <= high_limit:
                data_index[i] = 1
                keepers = keepers + 1
        if (data_in[i,4] - data_in[i,2]) < ug_limit:
            data_index[i] = 0
        #Remove stars randomly according to detection efficiency...
        app_mag = data_in[i,2] #(data_in[i,2] + 5.*(m.log10(new_dist*1000) - 1.))
        det_eff = s[0] / (np.exp(s[1]*(app_mag - s[2])) + 1.)
        check = np.random.uniform()
        #print app_mag, det_eff, check
        if (check > det_eff):
            data_index[i] = 0
            sig_boot = sig_boot + 1
    print '#-Stars in data set:', length, 'Stars in selection:', keepers
    print '#-Booted', sig_boot, 'stars due to detection efficiency'
    new_data = sc.zeros((sc.sum(data_index),width))
    j = 0
    for i in range(length):
        if data_index[i] == 1:
            new_data[j,:] = data_in[i,:]
            j = j + 1
    if (j != sc.sum(data_index)):  print '!!! stars out not equal to wanted #:', j, sc.sum(data_index)
    return new_data
    
    
def con_analysis(uncon_data, con_data, left_limit = 0.0, right_limit = 0.3):
    #unconvolved analysis
    l, w = uncon_data.shape
    o_right = 0
    o_left = 0
    o_middle = 0
    for i in range(l):
        if ( (uncon_data[i,2] - uncon_data[i,3]) < left_limit):
            o_left = o_left + 1
        elif ( (uncon_data[i,2] - uncon_data[i,3]) > right_limit):
            o_right = o_right + 1
        else:
            o_middle = o_middle + 1
    print "#-Before Convolution:  Left:  in region:  Right:"
    print o_left, o_middle, o_right
    #Convolved Analysis
    l, w = con_data.shape
    c_right = 0
    c_left = 0
    c_middle = 0
    for i in range(l):
        if ( (con_data[i,2] - con_data[i,3]) < left_limit):
            c_left = c_left + 1
        elif ( (con_data[i,2] - con_data[i,3]) > right_limit):
            c_right = c_right + 1
        else:
            c_middle = c_middle + 1
    print "#-After convolution:  Left:  in region:  Right:"
    print  c_left, c_middle, c_right
    #Finding differences
    delta_l = o_left - c_left      #These are fluxes, so they are reversed
    delta_m = c_middle - o_middle  #This is a standard "final-minus-initial"
    delta_r = o_right - c_right    #These are fluxes, so they are reversed
    print "#-Flux into f-turnoff region, from left:  Total change:  from right:"
    print delta_l, delta_m, delta_r
    return 1

""" Returns Poissonian error for a single bin count; either in simple form or in a rigorous
    approximation function, as provided in Gehrels 1986, for the 1-sigma confidence level
    NOTE:  Lower limit not yet implemented!  Can make much more rigorous!!!"""
def poisson_errors(n, rigorous=False):
    #upper error bar
    lambda_up = n + np.sqrt(n+1.0) + 1.0
    #lower error bar
    #lambda_down = n*( (1.0 - (1.0/(9.0*n)) - (1.0/(3.0*np.sqrt(n))))**3 )
    return lambda_up

''' Deprecated
def sigmoid_error(hist_in, real_dist, new_dist):  #CHANGE!
    s = [0.9402, 1.6171, 23.5877]
    hist_out = (sc.ones(hist_in.shape))*1000.0
    hist_out[:,0] = hist_in[:,0]
    for i in range(len(hist_in[:,0])):
        #Get abs mags of new 
        ori_mag = (hist_in[i,0] + 5.*(m.log10(real_dist*1000) - 1.))
        new_mag = (hist_in[i,0] + 5.*(m.log10(new_dist*1000) - 1.))
        #Calculate detection efficiency, as per Cole 2008
        ori_def = s[0] / (np.exp(s[1]*(ori_mag - s[2])) + 1.)
        new_def = s[0] / (np.exp(s[1]*(new_mag - s[2])) + 1.)
        #Claculate new, detection efficiency-corrected bin height
        hist_out[i,1] = (new_def / ori_def)*hist_in[i,1]
        print '#-bin mag', hist_in[i,0], 'to', new_mag,'modifier:', (new_def / ori_def)
    return hist_out
'''