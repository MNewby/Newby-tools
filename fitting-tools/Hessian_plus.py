#! /usr/bin/env python

import scipy as sc
import numpy as np
import subprocess as sp
import time

"""This program contains code for finding errors from R-squared functions.
Matthew Newby (RPI), Oct 16, 2010
http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node8.html
N_data = number of data points in original data
"""

class results():
    # Need function, best-fit data
    def __init__(self):
        self.params = []
        self.steps = []
        self.errors = []
	self.function = mw_func 

""" STUFF TO CHANGE FOR EACH RUN """
luafile = "ModfitParams.lua"
def mw_func(params):
    write_lua(params, luafile)
    sts1 = sp.call("./milkyway_separation -f -i -s stars-15.txt -a temp.lua", shell=True)
    likelihood = get_likelihood()
    #delete files
    sts2 = sp.call("rm temp.lua", shell=True)
    sts3 = sp.call("rm stderr.txt", shell=True)
    return likelihood

def get_hessian_errors(result, verbose=0, like=1):
    length = len(result.params)
    hessian = sc.zeros((length, length),float)
    base_params = sc.zeros(length, float)
    for i in range(length): base_params[i] = result.params[i]
    test_params = sc.zeros(length, float)
    #Build Hessian
    for i in range(length):
        for j in range(length):
            print "# - Finding element {0}, {1}".format(i, j)
            #H_1[i,j] = L(Q[j]+h[j], Q[i]+h[i])
            for k in range(length): test_params[k]=base_params[k]
            test_params[j] = test_params[j] + result.steps[j]
            test_params[i] = test_params[i] + result.steps[i]
            H_1 = result.function(test_params) #R_squared(function, test_params, x, y, sigma)
            if like==1:  H_1 = np.exp(-H_1/2.0)
            #H_2[i,j] = L(Q[j]-h[j], Q[i]+h[i])
            for k in range(length): test_params[k]=base_params[k]
            test_params[j] = test_params[j] - result.steps[j]
            test_params[i] = test_params[i] + result.steps[i]
            H_2 = result.function(test_params) #R_squared(function, test_params, x, y, sigma)
            if like==1:  H_2 = np.exp(-H_2/2.0)
            #H_3[i,j] = L(Q[j]+h[j], Q[i]-h[i])
            for k in range(length): test_params[k]=base_params[k]
            test_params[j] = test_params[j] + result.steps[j]
            test_params[i] = test_params[i] - result.steps[i]
            H_3 = result.function(test_params) #R_squared(function, test_params, x, y, sigma)
            if like==1:  H_3 = np.exp(-H_3/2.0)
            #H_4[i,j] = L(Q[j]-h[j], Q[i]-h[i])
            for k in range(length): test_params[k]=base_params[k]
            test_params[j] = test_params[j] - result.steps[j]
            test_params[i] = test_params[i] - result.steps[i]
            H_4 = result.function(test_params) #R_squared(function, test_params, x, y, sigma)
            if like==1:  H_4 = np.exp(-H_4/2.0)
            #H[i,j] = (H_1[i,j] - H_2[i,j]-H_3[i,j]+H_4[i,j]) / (4*h[i]*h[j])
            #h[k] is the step size for likelihood determination; use same as gradient?
            hessian[i,j] = (H_1 - H_2 - H_3 + H_4) / (4.0*result.steps[i]*result.steps[j])
    #check that it is symmetric
    for i in range(length):
        for j in range(length):
            if (hessian[i,j] != hessian[j,i]): 
                symmetric=False
                if verbose:  print '!!! Hessian not symmetric!'
            else:  symmetric=True
    #Convert to matrix type and invert
    hessian_matrix = np.matrix(hessian)
    print '#--Hessian Matrix:'
    print hessian_matrix
    hessian_inverse = hessian_matrix.I
    #read off diagonals and calculate errors
    errors = sc.zeros(length, float)
    for i in range(length):
        if (hessian_inverse[i,i] < 0.0):  print '!!!Hessian error', i, 'below zero!!!'
        errors[i] = np.sqrt(2.0*abs(hessian_inverse[i,i]))  #*(1.0/len(x))
    if verbose:  print '#---Hessian Errors:', errors
    return errors


def modify_steps(result, scale, tolerance=0.20):
    test=0
    for i in range(len(result.params)):
        diff = result.steps[i] - result.errors[i]
        if (abs(diff) >= (tolerance*result.params[i]) ):
            result.steps[i] = result.steps[i] - (diff*scale)
            test=test+1
    if test > 0:
        print "# - {0} parameters have not converged, starting next loop".format(test)
        return scale*0.8
    else:
        print "# - All parameters have converged.  Finishing..."
        return 0
# test &= (abs(result.steps[i] - result.errors[i]) < (tolerance*result.params[i]) ) 


def smart_Hessian(result, loops=10, scale=1.0, tolerance=0.2):
    t0 = time.time()
    while loops > 0:
        print "# - Starting Loop {0} (counting down); {1} seconds elapsed".format(loops, (time.time()-t0))
        result.errors = get_hessian_errors(result)
        scale = modify_steps(results, scale, tolerance)
        if scale==0:  break
        loops = loops-1
    if loops < 1:  print "!!! - Exited due to loop threshold"
    print "# - Hessian Errors:  {0}".format(result.errors)
    
    
def read_lua(filename):
    goodlines = [5,6, 11,12,13,14,15,16, 20,21,22,23,24,25, 29,30,31,32,33,34]
    params = []
    infile = open(filename, "r")
    number = 0
    for line in infile:
        number = number+1
        if number in goodlines:  params.append(float(line.split("=")[-1].strip().strip(',') ))
    infile.close()
    return params


def write_lua(params, inname, outname="temp.lua"):
    goodlines = [5,6, 11,12,13,14,15,16, 20,21,22,23,24,25, 29,30,31,32,33,34]
    commalines = [5, 11,12,13,14,15, 20,21,22,23,24, 29,30,31,32,33]
    infile = open(inname, "r")
    outfile = open(outname, "w")
    number, p = 0, 0
    for line in infile:
        number = number+1
        if number in goodlines:
            hold = line.split("=")
            if number in commalines:  out=str(params[p])+",\n"
            else:                     out=str(params[p])+"\n"
            outfile.write(hold[0]+" = "+out)
            p=p+1
        else:  outfile.write(line)
    infile.close()
    outfile.close()
    return 1

def get_likelihood(likefile="stderr.txt"):
    infile = open(likefile, "r")
    for line in infile:
        if line[0:19]=="<search_likelihood>":
            hold = line.split(" ")[1]
            break
    infile.close()
    return float(hold)

if __name__ == "__main__":
    result = results()
    result.params = read_lua(luafile)
    result.steps = [0.1,5.0, 0.2,20.0,10.0,1.0,1.0,1.0, 0.2,20.0,10.0,1.0,1.0,1.0, 0.2,20.0,10.0,1.0,1.0,1.0 ]
    results.errors = []
    results.function = mw_func
    smart_Hessian(result)
    #print get_likelihood()
    
