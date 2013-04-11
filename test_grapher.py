#! /usr/bin/env python

import scipy as sc
import math as m
import numpy as nu
import matplotlib as mpl
import matplotlib.pyplot as plt

"""This program is just a place to test new graphing functions.

Created by Matthew Newby (RPI), Dec 16, 2009.
"""

print 'test start successful'

"""This funcion generates a line with gaussian errors in the y coordinate.  
    
    Inputs:
        m, b: the parameters of a straight line: y=mx+b
        points: number of points to generate
        xstep: spacing between points in x-demension
        ysigma: standard deviation or y errors
        data: 1 if data file is to be saved, 0 otherwise
"""

def make_line(m=1.0, b=25.0, points=100, xstep=1.0, ysigma=5.0, data=1):
    
    #generate arrays
    x = sc.zeros(points)
    y = sc.zeros(points)
    y_center = sc.zeros(points)

    #initialize random seed
    nu.random.seed(10)
    
    #fill arrays
    for i in range(points):
        x[i] = i*xstep
    y_center = m*x + b
    y = (ysigma*nu.random.randn(points) ) + y_center
    print x,y
    
    #save array to file
    if data==1:
        data_out = raw_input('data file name, ending in .txt: ')
        if data_out == '':
            data_out = 'new_file.txt'
        for j in range(points):
            f = open(data_out, "w")
            f.write(str(x[j])+','+'\t'+str(y[j]))  
            #f.write(y[j])
        f.close()
    print 'file', data_out, 'successfully written and closed'

    #plot line
    plt.scatter(x, y)
    plt.show()

    print 'ending program...'
