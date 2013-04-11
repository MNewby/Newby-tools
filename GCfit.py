#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import gauss as ga

"""This program contains code for fitting gaussians to Globular cluster g histograms.

Matthew Newby (RPI), April 22, 2010
"""
#########################################################################################
GCname = "Pal 5 cluster"         	#Name of cluster being analyzed
myfile = "HR_Pal_5_cluster.csv"     	#Name of input .csv file
backfile = ".txt"
distance = 7.5                		#Distance to cluster in kpc
gcolumn = 0				#column number for g-magnitude values

#########################################################################################



    if Nhistgauss == 1:
        #x_hist, y_hist = normal_binned(col3, datarows)
        x_hist, y_hist = binned_background(col3, col7)
        nbins = len(x_hist)
        fit_params = ga.gauss_fit(x_hist, y_hist, nbins, 1,1)
        axis = 0.0 #(np.ma.max(x_hist) - np.ma.min(x_hist))/2.0 #non-zero lengthens x-axis
        x_line = sc.arange( (np.ma.min(x_hist)-axis), (np.ma.max(x_hist)+axis), 0.01)
        y_line = sc.zeros(len(x_line))
        mu_apparent = fit_params[0]
        fit_params[0] = fit_params[0] - 5.*(m.log10(distance*1000) - 1.)
        mu_str = '$\mu_M=$ ' + str(fit_params[0])
        sigmal_str = '$\sigma_l=$ ' + str(fit_params[1])
        sigmar_str = '$\sigma_r=$ ' + str(fit_params[2])
        A_str = '$A=$ ' + str(fit_params[3])
        mug_str = '$\mu_g=$ ' + str(mu_apparent)
        print '-best fit parameters:', fit_params
        print '-apparent magnitude mu:', mu_apparent
        fwhm = ga.full_width_half_max(fit_params[1], fit_params[2])
        FWHM_str = '$fwhm=$ ' + str(fwhm)
        print '-full width at half max: ', fwhm
        #Shifts x values to absolute magnitudes, creates a gaussian profile from params
        for j in range(len(x_hist)):
            x_hist[j] = x_hist[j] - 5.*(m.log10(distance*1000) - 1.)
        for j in range(len(x_line)):
            x_line[j] = x_line[j] - 5.*(m.log10(distance*1000) - 1.)
            if x_line[j] < fit_params[0]:
                sigma = fit_params[1]
            else:
                sigma = fit_params[2]
            exponent = -0.5 * ( (x_line[j] - fit_params[0]) / sigma )**2
            stuff = m.exp(exponent)
            y_line[j] = stuff*fit_params[3]
        #Creates a line to simulate binning
        x_bins = sc.arange( (np.ma.min(x_hist)-axis), (np.ma.max(x_hist)+axis), 0.01)
        y_bins = sc.zeros(len(x_bins))
        span = x_hist[1] - x_hist[0]
        bin_min = x_hist[0] + (span/2.0) #Starts at first bin-bin changover
        k=0
        for j in range(len(x_bins)):
            if (x_bins[j] > (bin_min + (span*k))):
                k = k + 1
            y_bins[j] = y_hist[k]

        #Shifts original x values to absolute magnitudes, prepares them for binning
        #x_shift = sc.zeros(datarows)
        #x_shift = col3 - 5.0*(m.log10(distance*1000) - 1.)
        plt.figure(1)
        plt.plot(x_bins, y_bins, 'b-')
        plt.scatter(x_hist, y_hist, marker='d' )
        #plt.hist(x_shift, nbins, histtype='step') 
        plt.plot(x_line, y_line, 'g-')
        plt.title('line fit of ' + GCname)
        plt.xlabel('Absolute Magnitude')
        plt.ylabel('Counts')
        if (plot_params == 1):
            txt_string = mu_str+'\n'+sigmal_str+'\n'+sigmar_str+'\n'+A_str+'\n'+mug_str+'\n'+FWHM_str
            plt.text(np.ma.min(x_line), (np.ma.max(y_line)*0.6), (txt_string))
        plt.show()

