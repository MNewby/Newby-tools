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

'''python script for running gc scripts
Matthew Newby, Sept 6, 2010'''

'''Add Limits to parameters and data! '''

"""
# Generate Gaussian Data
sig = 10.0
x = sc.arange(0.0, 10.0, 0.1)
y = (func.gaussian_function(x, [10.0, 3.5, 0.8]) )#+ np.random.normal(0.0, sig, len(x)))
out = sc.zeros((len(x),3))
for i in range(len(x)):
    out[i,0], out[i,1], out[i,2] = x[i], y[i], sig
f.write_data(out, "test_out_p.txt")
"""

#MAIN
#Setup variables
filename = 'Pal13_backsub_cluscut.txt'
function = func.get_2gauss_y
identifier = ['text', 'A', 'B', 'C', 'D', 'F', 'G', 'H']

#Load data
suffix = filename[-4:]
if suffix == '.txt':
    data = f.read_data(filename)
elif suffix == 'text':
    data = f.read_data(filename, ',')
elif suffix == '.csv':
    data = f.read_data(filename, ',')    
else:
    print "WHAT?  I don't know how to read", suffix
#data = make_line_errors(10.0, 20.0, yerr=15.0)
#data = make_gauss_errors(2.0, 2.0, yerr=0.0)

#Set data
#columns in data file for respective values
x_col, y_col, sig_col = data[:,1], data[:,0], data[:,1] #sc.ones(len(data[:,0]),float) #
#y_col, edges = sc.histogram(data[:,15], bins=20, range=(-200.0, 200.0))
#x_col = (edges + 10.0)[:-1] #sc.arange( (low+(size/2.0)), 200.0, size)
sig_col = func.poisson_errors(y_col)

""" Initializations - choose correct number of parameters"""
""" 6 params """
#start_params = sc.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  #A function that gives a good start is nice here
#steps = sc.array([0.01,0.01,0.01,0.01,0.01,0.01])  #1/10 of expected parameter order is usually good
""" 4 params """
start_params = sc.array([4.456, 0.814, 0.226, 8.43])  #A function that gives a good start is nice here
steps = sc.array([0.001,0.001,0.001,0.001])  #1/10 of expected parameter order is usually good
""" 3 params """
#start_params = sc.array([5.0, 5.0, 1.5])  #A function that gives a good start is nice here
#steps = sc.array([0.001,0.001,0.001])  #1/10 of expected parameter order is usually good
""" 2 params """
#start_params = sc.array([1.0, 0.001])  #A function that gives a good start is nice here
#steps = sc.array([0.01,0.0001])  #1/10 of expected parameter order is usually good
""" 1 params """
#start_params = sc.array([60.0])  #A function that gives a good start is nice here
#steps = sc.array([0.01])  #1/10 of expected parameter order is usually good

""" fit data """
MCMC_fit, errors, best_fit = mc.do_MCMC(function, start_params, steps, x_col, y_col,
                                                                      sig_col, name = identifier, number_steps=10000, save=1)
best_fit = gd.gradient_descent(function, best_fit, x_col, y_col, sig_col)
#best_fit = gd.gradient_descent(function, start_params, x_col, y_col, sig_col, to_screen=1)
print '#---GD Best fit:', best_fit

#plot data and fit
header = [(identifier[0]+str(best_fit)), identifier[1], identifier[2]]
if (pdf.plot_function(x_col, y_col, function, best_fit, save_name=identifier[0], 
                      title=[header[0], r'$radius$', r'$counts$'], save=1) == 1):
    print '#---Plotting successful'

#Get Errors
errors = he.get_hessian_errors(function, best_fit, x_col, y_col, sig_col, (steps) )

#y_err=sig_col,
        
print '#---All tasks complete'
