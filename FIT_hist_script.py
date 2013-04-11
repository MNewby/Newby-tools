import math as m
import numpy as np
import scipy as sc
import files as f
import gradient_descent as gd
import monte_carlo as mc
import Hessian_errors as he
import plot_data_function as pdf
import time
import functions as func
import glob

'''python script for running gc scripts
Matthew Newby, Sept 6, 2010'''

'''Add Limits to parameters and data! '''

#MAIN
#Setup variables
filename = glob.glob("./An_Project/*real_noG1_20.0*.txt")
print filename
function = func.double_gaussian

for i in range(len(filename)):
    #Load data
    uid = "fit_"+filename[i][-13:-4]
    identifier = [uid, 'A', 'B', 'C', 'D', 'F', 'G', 'H']
    data = f.read_data(filename[i])

    #Set data
    #columns in data file for respective values
    x_col, y_col, sig_col = (data[:,1]+0.05), data[:,0], func.poisson_errors(data[:,0])
 
    """ Initializations - choose correct number of parameters"""
    """ 6 params """
    start_params = sc.array([5.0, -1.0, 1.0, 2.5, 1.0, 0.5])  #A function that gives a good start is nice here
    steps = sc.array([0.1,0.01,0.01,0.1,0.01,0.01])  #1/10 of expected parameter order is usually good

    """ fit data """
    MCMC_fit, errors, best_fit = mc.do_MCMC(function, start_params, steps, x_col, y_col,
                              sig_col, name = identifier, number_steps=10000, save=1)
    best_fit = gd.gradient_descent(function, best_fit, x_col, y_col, sig_col)
    print '#---GD Best fit:', best_fit

    #plot data and fit
    header = [(identifier[0]+str(best_fit)), identifier[1], identifier[2]]
    if (pdf.plot_function(x_col, y_col, function, best_fit, save_name=identifier[0], 
                          title=[header[0], r'$radius$', r'$counts$'], save=1) == 1):
        print '#---Plotting successful'

    #Get Errors
    errors = he.get_hessian_errors(function, best_fit, x_col, y_col, sig_col, (steps*0.1) )
    print "$ {0}".format(uid)
    print "$ {0}".format(best_fit)
    print "$ {0}".format(errors)
    
print '#---All tasks complete'
