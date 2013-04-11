#! /usr/bin/env python  #'Bang' line - modify as needed

import sys
import math as m
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')  #Use this to save plots as postscripts; not sure if necessary
import matplotlib.pyplot as plt

"""This is a sample script for plotting data using matplotlib, which requires
scipy and numpy to run.

Matthew Newby (RPI), Jan 5, 2010
"""

def make_data():
    #Set parameters for sine wave functions:
    A_1, A_2 = 1.0, 1.0  #Wave Amplitudes
    T_1, T_2 = 1.0, 2.5  #Wave Periods
    phi_1, phi_2 = 0.0, 0.0  #Phase values
    k_1, k_2 = ( (2.0*np.pi) / T_1), ((2.0*np.pi) / T_2)  #Clean up the coeffiecients
    #Create sine waves and add them:
    t = np.arange(0.0,10.0,0.01) #sets t to an array of values from 0.0 to 10.0, in increments of 0.01
    x_1 = A_1*np.sin(k_1*t + phi_1)  #since 'sin()' is a numpy ('np') function, element-wise evaluation...
    x_2 = A_2*np.sin(k_2*t + phi_2)  #...is implied, so x_1 and X_2 are automatically arrays.
    return t, x_1, x_2  #returns a tuple of t and sine wave values along t

def plot():  #Change so that your arguments are passed in correctly
    """Plot Initializations"""
    plt.figure(1)
    plt.title('Sine Wave Addition', fontsize=10)
    plt.xlabel('$t$') #The dollar signs produce latex-formatted text (looks pretty for variables)
    plt.ylabel('$x$')
    """Get the Data"""
    x, y_1, y_2 = make_data()
    """Now to create our line plots"""
    plt.plot(x, (y_1+y_2), 'k-')  #Creates a solid black line for the added waves
    plt.plot(x, y_1, 'g:')  #Creates a green dotted line for sine wave 1
    plt.plot(x, y_2, 'b:')  #Creates a blue dotted line for sine wave 1
    """Let's look at the plot"""
    plt.show()
    return 1

if __name__ == '__main__':
    plot()
    print "Done Plotting"