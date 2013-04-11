#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_ex():

    plt.figure(1)  		#create figure
    plt.subplot(221)		#create subplots, rows, columns, plot number
    plt.scatter(x, y)		#plotting function, 
    plt.axis('normal')          #axis options; scaling, etc.
    plt.title('title')   	#plot title
    plt.xlabel('x_axis')	#x axis title
    plt.ylabel('y_axis')	#y_axis title
    plt.xlim(x_min, x_max)	#x axis range
    plt.ylim(y_min, y_max)	#y axis range
    plt.show()			#displays plots

    return 0
