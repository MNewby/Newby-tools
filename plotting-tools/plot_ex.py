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

Matthew Newby (RPI), Dec 10, 2010
"""

def plot(argv):  #Change so that your arguments are passed in correctly
    """Plot Initializations - can also create subplots if you want to be fancy -
        see online documentation"""
    plt.figure(1)
    
    """Text Values - note that these accept latex-formatted strings!  Also, any
        command that produces text can take 'text' arguments:
        http://matplotlib.sourceforge.net/api/artist_api.html#matplotlib.text.Text"""
    plt.title('I am a plot!')
    plt.xlabel('I am $x$!')
    plt.ylabel('I am $y$!')
    
    """Plotting Calls, see http://matplotlib.sourceforge.net/index.html for a
        list of functions, and individual functions for all possible arguments"""
    plt.plot(x, y)  #creates a line plot
    plt.scatter(x, y)  #creates a scatter plot
    plt.errorbar(x, y, xerr=None, yerr=None)  #creates a line plot with error bars
    #- can be modifed to act like 'plot' or 'scatter', so this guy is really all you need...
    plt.hist(x, bins=10)  #makes a histgram of x.
    #I recommend setting your bin size then finding the number of bins from there...
    
    """Plot limits setting commands - if not called, matplotlib will set the limits to
        fit all plotted values, and scale the axes to make the final plot square"""
    plt.xlim(x_min, x_max)  
    plt.ylim(y_min, y_max)
    plt.axis('scaled')  #Scales the plot to make the axes proportional

    """Sets the tick mark locations and labels, inputs are arrays.  If not called, then
        matplotlib will choose what it thinks is best, Which is usually not bad.  Also, if
        you call these functions without arguments, they return a tuple of current
        'locations' and 'labels' values"""
    plt.xticks(locations, labels)  
    plt.yticks(locations, labels)  

    """Use this if you want to have your code saved to a file"""
    plt.savefig(save_name, papertype='letter')  #saves as an 8.5"x11" postscript

    """Use this to just show the plot using the helpful GUI - do only one of either
        this or savefig, though!""" 
    plt.show()

    """It is always a good idea to clean up after yourself, especially if you still
        want to make more plots"""
    plt.close('all')
    
    """Let's get the hell outta here!"""
    return 1


if __name__ == '__main__':
    plot()