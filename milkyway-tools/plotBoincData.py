#! /usr/bin/env python  #'Bang' line - modify as needed

import math as m
import numpy as np
import scipy as sc
import sys
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt

"""Plots data from a series of BOINC stream-integration runs.  Input data must be
separated by spaces, but can be broken up onto different lines in any fashion.
Inputs:
1 - File to be read
2 - number of streams in data
3 - (optional) unique identifier for output files
Matthew Newby (RPI), Jan 20, 2011
"""

# FROM "files.py"
"""Builds a list, where each member is the data from a single run.
Each data set, for n streams:
Likelihood[0], q[1], R0[2], epsilon_1[3], mu_1[4], R_1[5], theta_1[6], phi_1[7], sigma_1[8], ...,
epsilon_n[3+(6*n)], mu_n[4+(6*n)], R_n[5+(6*n)], theta_n[6+(6*n)], phi_n[7+(6*n)], sigma_n[8+(6*n)]"""
def read_boinc_results(filename, num_streams):
	readfile = open(filename, "r")
	#output data list, run list, iteration counter, # parameters per run
	data, run, i, parameters = [], [], 1, (3+(6*num_streams))
	for line in readfile:
		if i > parameters:
			data.append(run)
			run = []
			i = 1
		if (line.strip() == ''): continue
		if (line[0] == "#"): continue
		holder = line.split()
		for param in holder:
			run.append(param.strip(' ,'))
			i=i+1
	data.append(run)
	readfile.close()
	print '#---Boinc separation file', filename, 'with', num_streams, 'streams successfully read.'
	return sc.array(data)

def plot_boinc(args):
    if len(args) < 3:  identifier = 'new'
    else: identifier = str(args[2])
    data_in = read_boinc_results(args[0], int(args[1]))
    figs = []
    x = [0, 10, 20, 30, 40]
    Names = [['Fitness', '$q$', '$R_0$'], ['$\epsilon$', '$\mu$',
                'radius, $R$', '$\theta$', '$\phi$', '$\sigma$']]
    for i in range(int(args[1])+1):
        fig = plt.figure(i)
        if i == 0:  #Fitness, background plots
            plt.title('Background')
            for j in range(3):
                plt.subplot(3,1,(j+1))
                plt.errorbar(x, data_in[:,j], fmt='b-', color='black', marker ='s', ms=4)
                plt.xticks(fontsize=8)
                plt.yticks(fontsize=8)
                plt.ylabel(Names[0][j], fontsize=6)
            save_name = 'Background_'+identifier+'.ps'
        else:
            plt.title( ('Stream '+str(i) ))
            for j in range(6):
                plt.subplot(6,1,(j+1))
                plt.errorbar(x, data_in[:,(3+j+(6*(i-1)))], fmt='b-', color='black', marker ='s', ms=4)
                plt.xticks(fontsize=6)
                plt.yticks(fontsize=6)
                plt.ylabel(Names[1][j], fontsize=6)
            save_name = 'Stream'+str(i)+'_'+identifier+'.ps'
        figs.append(fig)
        plt.savefig(save_name, papertype='letter')
        print "#---", save_name, "plotted and saved"
    return 0

if __name__ == "__main__":
    plot_boinc(sys.argv[1:])