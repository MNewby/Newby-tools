import csv
import math as m
import numpy as np
import scipy as sc
import files as f
import time
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt

"""This script fits a line to a subset of isocrone data.
Extra crap is at the bottom, in case an actual maximum likelihood
fit is needed.

  Just plots data with fit function right now.
Matthew Newby, RPI, Aug 26, 2010"""

def Gradient_Descent(data_in, st_parameters=[0.0,0.0], iterations=1000, precision=0.000000001):
    #Modifying this to fit isocrones to fiducial sequences.
    #May produce bad results if Amp. is small; finds strange local min.
    #Params are [A,B] for y = Ax+b
    #for manual starting point, beware of new params being better (use new_params instead):
    print '-Starting Gradient Descent...'
    print st_parameters
    old_params = [100.0, 100.0]
    new_params = st_parameters
    gradient = [1000.0, 1000.0]
    RR_old = r_squared_linear(data_in[:,0], (data_in[:,2]-data_in[:,1]), \
                              data_in[:,3], old_params)
    RR_new = r_squared_linear(data_in[:,0], (data_in[:,2]-data_in[:,1]), \
                              data_in[:,3], new_params)
    loop = 0
    overrun = 0
    scale = 0.00000001
    print RR_old, RR_new, precision
    while (abs(RR_old-RR_new) > precision):  #MAY STILL MISS LAST STEP!???
        if (RR_new < RR_old):
            print '### I moved!  r-squared change:', abs(RR_old - RR_new)
            RR_old = RR_new
            for k in range(len(old_params)):
                old_params[k] = new_params[k]
            scale = scale*1.5
        else:
            scale = scale*0.5
        loop = loop + 1
        if (loop > iterations):
            print '-Exited by passing iteration threshold'
            overrun = 1
            break
        gradient[0] = dR2_dA(data_in[:,0], (data_in[:,2]-data_in[:,1]), data_in[:,3], old_params)
        gradient[1] = dR2_dB(data_in[:,0], (data_in[:,2]-data_in[:,1]), data_in[:,3], old_params)
        print 'gradient', gradient
        for i in range(len(old_params)):
            new_params[i] = old_params[i] - scale*gradient[i]
        RR_new = r_squared_linear(data_in[:,0], (data_in[:,2]-data_in[:,1]), \
                              data_in[:,3], new_params)
        print 'scale', scale, 'params', old_params, loop
    if overrun == 0:
        print '-Gradient Descent successful, exited after ', loop, ' iterations'
    return old_params

def dR2_dA(x, y, sigma, params):
    A_gradient = 0.0
    for i in range(len(x)):
        peice = -1.0*x[i]*(y[i]-(params[0]*x[i])-params[1])
        #peice = peice / (sigma[i]*sigma[i])
        A_gradient = A_gradient + peice
        #print 'delA', A_gradient
    return (2*A_gradient)

def dR2_dB(x, y, sigma, params):
    B_gradient = 0.0
    for i in range(len(x)):
        peice = -1.0*(y[i]-(params[0]*x[i])-params[1])
        #peice = peice / (sigma[i]*sigma[i])
        B_gradient = B_gradient + peice
        #print 'delB', B_gradient
    return (2*B_gradient)

def r_squared_linear(x, y, sigma, params):
    summation = 0.0
#    print x, y, sigma
#    for i in range(len(x)):
    part = ( (y - (x*params[0]) - params[1]) )**2  #deleted sigma for testing purposes
    #print part
    summation = sc.sum(part)
    #print 'sum', summation
    #sigma might be bad?  should always be greater than or equal to 1.0?
    return summation


"""MAIN"""
in_file = 'NGC_7078_KI03_out.txt'
A = -0.00648959
B = 0.05001323
limits = [0.0,8.0]
#remember - 'y' is [delta]g-r [2-1], 'x' is Mg [0], sigmas are [3]
#fitting 'y = Ax+B'
in_data = f.read_data(in_file)
fit_l, fit_w = in_data.shape
#Get only the wanted data
keepers = 0
for i in range(fit_l):
    if ( (in_data[i,0] > limits[0]) and (in_data[i,0] < limits[1])):
        keepers = keepers + 1
fit_data = sc.zeros((keepers, fit_w),float)
j = 0
for i in range(fit_l):
    if ( (in_data[i,0] > limits[0]) and (in_data[i,0] < limits[1])):
        fit_data[j,:] = in_data[i,:]
        j = j + 1
fit_l, fit_w = fit_data.shape
"""#Starting Params
np.random.seed(int(time.time() ) )
start_params = [0.0,0.0]
start_params[0], start_params[1] = np.random.random_integers(-100,100,2)
#Fit with gradient descent
final_params = Gradient_Descent(fit_data, start_params)
A = final_params[0]
B = final_params[1]"""
plt.figure()
#plt.scatter((fit_data[:,2] - fit_data[:,1] +(A*fit_data[:,0])+B), fit_data[:,0], \
#    c='black', marker='+', label='subtracted data')
#plt.plot([0.0,0.0], [np.ma.min(fit_data[:,0]),np.ma.max(fit_data[:,0])], 'r-', label='_nolegend_')
#plt.errorbar(fit_data[:,2],fit_data[:,0],xerr=fit_data[:,3], fmt=None, ecolor='g')
plt.scatter(fit_data[:,2],fit_data[:,0], c='green', label='fiducial sequence')
plt.scatter((fit_data[:,1]+(A*fit_data[:,0])+B),fit_data[:,0], label='Girardi isocrone')
plt.scatter((fit_data[:,2] - fit_data[:,1]), fit_data[:,0], marker='+')
#plt.scatter(in_data[:,0], (in_data[:,2] - in_data[:,1]), c='g')
plt.plot([((A*limits[0])+B), ((A*limits[1])+B)], (limits))
plt.ylim(limits[1], limits[0])
plt.show()