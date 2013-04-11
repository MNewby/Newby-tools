#! /usr/bin/env python

import files as f
import hist as h
import scipy as sc
import math as m
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import functions as func
import gradient_descent as gd
import monte_carlo as mc
import Hessian_errors as he
import plot_data_function as pdf

"""This script is for simple scripts.

Matthew Newby (RPI), May 4, 2011
"""

spread = 0.5

skip = [0,0,0,0,0,0,0,0,0,0]
datafile = ['noU_NGC_4147_cluster.csv',
            'noU_NGC_5024_cluster.csv',
            'noU_NGC_5053_cluster.csv',
            'noU_NGC_5272_cluster.csv',
            'noU_NGC_5466_cluster.csv',
            'noU_NGC_5904_cluster.csv',
            'noU_NGC_6205_cluster.csv',
            'noU_NGC_7078_cluster.csv',
            'noU_NGC_7089_cluster.csv',
            'noU_Pal5_cluster.csv',
           ]
DISTANCE = [19.3, 18.7, 18.5, 10.4, 15.6, 8.0, 7.7, 11.0, 11.5, 21.0]

results_file = 'results_new_corrected_singles.txt'

results = f.read_data(results_file)

turnoff = []
for i in range(len(datafile)):
    if skip[i] == 1:  continue
    data = f.read_data(datafile[i], ',')
    mu = results[i,2] + 5.*(m.log10(DISTANCE[i]*1000) - 1.)
    holder=[]
    for j in range(len(data[:,0])):
        if data[j,2] > (mu+spread):  continue
        if data[j,2] < (mu-spread):  continue
        holder.append(data[j,2]-data[j,3])
    histogram = h.make_hist(holder, 0.01)
    h.plot_histogram(histogram, 0.01, name=('hist_'+str(i)), x_label=r'$(g-r)_0$')
    #h.plot_multiple_hist(histogram[:,0], [histogram[:,1]], 0.01, name=('hist_'+str(i)), x_label=r'$(g-r)_0$')
    #x_col, y_col, sig_col = histogram[:,0], (histogram[:,1]/len(holder)), (func.poisson_errors(histogram[:,1])/len(holder) )
    x_col, y_col, sig_col = histogram[:,0], (histogram[:,1]), (func.poisson_errors(histogram[:,1]))
    #print x_col, y_col, sig_col
    start_params = sc.array([0.2, 0.1, sc.ma.max(y_col)])
    steps = sc.array([0.001, 0.001, 0.01])
    function = func.gaussian_function
    #Fit function
    MCMC_fit, errors, best_fit = mc.do_MCMC(function, start_params, steps, x_col, y_col,
                              sig_col, name=('histfit_'+str(i)), number_steps=10000, save=0)
    best_fit = gd.gradient_descent(function, best_fit, x_col, y_col, sig_col)
    errors = he.get_hessian_errors(function, best_fit, x_col, y_col, sig_col, steps)
    turnoff.append([best_fit[0], best_fit[1], errors[0], errors[1]])
#    if (pdf.plot_function(x_col, y_col, function, best_fit, save_name=('histfitplot_'+str(i)), 
#                      title=['', r'$(g-r)_0$', 'Normalized Counts'], save=1) == 1):
    plt.figure()
    plt.scatter(x_col, y_col)
    func_x = sc.arange( (sc.ma.min(x_col)*0.9), (sc.ma.max(x_col)*1.1), 0.01)
    func_y = function(func_x, best_fit)
    plt.plot(func_x, func_y)
    plt.savefig(('histfitplot_'+str(i)+'.ps'), papertype='letter')
    print  "#-Fit dataplot saved"
#print turnoff
print turnoff
if f.write_data(sc.array(turnoff), 'turnoff_out_fixed.txt') == 1:
    print "#-All Results Successfully saved"