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
# nohup python Hessian_plus.py &> dump.out &

sts0 = sp.call("rm temp.lua", shell=True)
sts00 = sp.call("rm stderr.txt", shell=True)

class results():
    # Need function, best-fit data
    def __init__(self):
        self.params = []
        self.steps = []
        self.errors = []
        self.function = mw_func
    def speak(self):
        #printa("Params: "+str(self.params)) 
        printa("Steps: "+str(self.steps))
        printa("Errors: "+str(self.errors))

""" STUFF TO CHANGE FOR EACH RUN """
dumpfile = "dump.out"
sts000 = sp.call("rm {0}".format(dumpfile), shell=True)
luafile = "params15.lua"
def mw_func(params):
    write_lua(params, luafile)
    sts1 = sp.call("./milkyway_separation -f -i -s stars-15.txt -a temp.lua", shell=True)
    likelihood = get_likelihood()
    #delete files
    sts2 = sp.call("rm temp.lua", shell=True)
    sts3 = sp.call("rm stderr.txt", shell=True)
    return likelihood

def do_MCMC(result, n_steps=1000):
    #function, init_params, step_sizes, x, y, sigma, name, number_steps=1000, save=1):
    #name is a list of strings of this form: run name, param name 1, param name 2, ...
    np.random.seed( int(time.time()) )
    lp, moves, positions = len(result.params), 0, []  #lp=length of parameters
    best_RR = 10000000.0 
    best_params = sc.zeros(lp)
    current_params = sc.zeros(lp)
    for i in range(lp):  current_params[i] = result.params[i]
    new_params = sc.zeros(lp)
    current_RR = result.function(current_params) #R_squared(function, current_params, x, y, sigma)
    if current_RR < 0.0:  current_RR = -1.0*current_RR  #fix for likelihoods
    for i in range(n_steps):
        #Record position
        positions.append(current_params.tolist())
        #Take a step
        for j in range(lp):
            new_params[j] = np.random.normal(current_params[j], result.steps[j])
        #Decide whether to move or not
        new_RR = result.function(new_params) #R_squared(function, new_params, x, y, sigma)
        if new_RR < 0.0:  new_RR = -1.0*new_RR
        compare = (current_RR / new_RR) 
        #if (np.random.uniform() < (np.exp(-1.0*new_RR) / np.exp(-1.0*current_RR) ) ):
        if (np.random.uniform() < compare):
            moves = moves + 1
            current_RR = new_RR
            for j in range(lp):
                current_params[j] = new_params[j]
        printa("{0} : {1}".format(current_params, current_RR))
        #record best fitness
        if (new_RR < best_RR):
            best_RR = new_RR
            for j in range(lp):
                best_params[j] = new_params[j]
        #if (i % 1000) == 0:
        #    print 'At iteration', i, current_params
    #make histogram 
    #centers, deviations = plot_MCMC_hist(positions, name, save)
    #print '#---Mean parameters:', centers, 'Parameter deviations:', deviations
    #print '#---Best parameters:', best_params, best_RR
    #print '#---After:', moves, ' moves out of ', number_steps,'iterations'
    #plot fit with function --- with both best and means?
    #if (plot_function(x, y, function, params, x_err=None, y_err=sigma,
    #              title=name[0], save=0, save_file=(name[0]+'_func_plot.ps') ) ==1):
    #    print '#---data and MCMC function fit successfully plotted'
    pos_array = sc.array(positions)
    printa("Original Parameters: {0}".format(result.params))
    printa("New Best?: {0} : {1}".format(best_params, best_RR))
    for i in range(lp):
        printa("Parameter {0}: {1}, {2}".format(i,sc.mean(pos_array[100:,i]),
                sc.std(pos_array[100:,i]) ) )
    return -1 #centers, deviations, best_params


def get_hessian_errors(result, verbose=1, like=1):
    length = len(result.params)
    hessian = sc.zeros((length, length),float)
    base_params = sc.zeros(length, float)
    for i in range(length): base_params[i] = result.params[i]
    test_params = sc.zeros(length, float)
    #Build Hessian
    for i in range(length):
        for j in range(length):
            #`print "# - Finding element {0}, {1}".format(i, j)
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
                if verbose:  printa('!!! Hessian not symmetric!')
            else:  symmetric=True
    #Convert to matrix type and invert
    hessian_matrix = np.matrix(hessian)
    printa('#--Hessian Matrix: \n{0}'.format(hessian_matrix))
    hessian_inverse = hessian_matrix.I
    #read off diagonals and calculate errors
    errors = sc.zeros(length, float)
    for i in range(length):
        if (hessian_inverse[i,i] < 0.0):  printa('!!!Hessian error '+str(i)+' below zero!!!')
        errors[i] = np.sqrt(2.0*abs(hessian_inverse[i,i]))  #*(1.0/len(x))
    if verbose:  printa('#---Hessian Errors: {0}'.format(errors))
    return errors


def modify_steps(result, scale, drag=0.9, tolerance=0.20):
    test=0
    for i in range(len(result.params)):
        diff = result.steps[i] - result.errors[i]
        if (abs(diff) >= (tolerance) ):
            result.steps[i] = result.steps[i]*(1.0 - scale*(diff/abs(diff)))
            #result.steps[i] = result.steps[i] - (diff*scale)  # Old, shitty way
            test=test+1
    if test > 0:
        printa("# - {0} parameters have not converged, starting next loop".format(test))
        return scale*drag
    else:
        printa("# - All parameters have converged.  Finishing...")
        return 0
# test &= (abs(result.steps[i] - result.errors[i]) < (tolerance*result.params[i]) ) 


def smart_Hessian(result, loops=10, scale=0.5, tolerance=20.0):
    """ This method is actually quite dumb """
    t0 = time.time()
    while loops > 0:
        printa("# - Starting Loop {0} (counting down); {1} seconds elapsed".format(loops, (time.time()-t0)))
        result.speak()
        result.errors = get_hessian_errors(result)
        scale = modify_steps(result, scale, tolerance)
        #result.speak()
        if scale==0:  break
        loops = loops-1
    if loops < 1:  printa ("!!! - Exited due to loop threshold")
    printa("# - Hessian Errors:  {0}".format(result.errors))
    printa("# - Total Time elapsed: {0} seconds".format(time.time()-t0))
    
    
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
    hold = -99999.9
    for line in infile:
        if line[0:19]=="<search_likelihood>":
            hold = line.split(" ")[1]
            break
    infile.close()
    return float(hold)

def printa(string, filename=dumpfile):
    outfile = open(filename, 'a')
    outfile.write(string+'\n')
    outfile.close()

if __name__ == "__main__":
    t0 = time.time()
    result = results()
    result.params = read_lua(luafile)
    result.steps = [0.1,1.0, 
        0.1,10.0,1.0,0.2,0.2,0.5, 
        0.1,10.0,1.0,0.2,0.2,0.5, 
        0.1,10.0,1.0,0.2,0.2,0.5 ]
    results.errors = []
    results.function = mw_func
    #smart_Hessian(result, loops=10) #change this line
    #print get_likelihood()
    do_MCMC(result, n_steps=1000)
    printa("Minutes Elapsed: {0}".format((time.time()-t0)/60.0) )
