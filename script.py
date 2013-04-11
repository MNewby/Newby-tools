import math as ma
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt

import plot_pack as pp

data2 = np.loadtxt("NGC_6205.csv", delimiter=",")

x = data2[:,2]-data2[:,3]
y = data2[:,2]

gerr = 3.11E-4 + np.exp(0.79*data2[:,2] - 20.0)
rerr = -2.63E-5 + np.exp(0.80*data2[:,3] - 19.8)
xerr = np.sqrt(gerr*gerr + rerr*rerr)
yerr = gerr


#Initializations - x and y are lists (or arrays) of data x,y points, sizes are bin sizes,
#    areas are (min, max) to be binned, None gives full range of data.
tester = pp.HistMaker(x, y, xsize=0.01, ysize=0.1, xarea=None, yarea=(15.0, 25.0) )
tester.cmap = 'color'             # Sets colormap;  'color' and 'bw' are hard-coded for convienience.
tester.varea = (0.0,None)         # Sets contrast range for pixels
tester.scale = None               # Special scaling, such as 'sqrt' or 'log'
tester.labels = [r'$(g-r)_0$',r'$g$']      # Labels for x,y axes
tester.xflip = 0   # Set=1 to flip x axis
tester.yflip = 1   # Set=1 to flip y axis

#tester.gblur(blur=1)  # Uncomment this to apply a Gaussian blur to the image

# Use this next line to apply an error convolution to the bins; requires x,y errors for each point.
#tester.H = pp.ErrBlurPoint(x, y, tester, con=30, xerror=xerr, yerror=yerr)

tester.plot()
