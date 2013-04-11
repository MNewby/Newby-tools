#! /usr/bin/env python

import files as f
import scipy as sc
import math as m
import numpy as np
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt


"""This program plots data sets of isochrones with fiducial sequences,
looking at the turnoff only.  OLD, not used.

Matthew Newby (RPI), Oct 2, 2010
"""

save = 1
#Could roll things up in lists.

def gaussian_function (x_in, params):
    exponent = -0.5*( (x_in - params[0]) / params[1])**2
    factor = 1.0 /(np.sqrt(2.0*m.pi*params[1]*params[1]))
    y_out = factor*np.exp(exponent)
    return y_out

"""Get Data"""
save_name = 'iso_compare_6205_clean3.ps'
data = f.read_data('NGC_6205_CG97_out.txt')

"""Setup Data for Plotting"""
y_limits = [7.5, 3.5]
iso_x = data[:,1]
y = data[:,0]
fid_x = data[:,2]
fid_err = data[:,3]
A = -0.01463
B = 0.08929
x = A*y + B

"""Plotting Functions"""
plt.figure()
plt.scatter((fid_x[:]-iso_x[:]),y, c='black', marker='+', label='subtracted data')
plt.plot([0.0,0.0], [y_limits[0],y_limits[1]], 'r-', label='_nolegend_')
plt.plot(x, y, 'b:', label='Line Fit')
plt.errorbar(fid_x[:],y,xerr=fid_err[:], fmt='o', c='b', ecolor='b', label='Fiducial Sequence')
plt.scatter(iso_x[:], y, c='r', marker='d', label='Girardi Isocrone')
plt.plot((fid_x - x), y, 'b:')
#plt.title(plot_title, fontsize='small')
plt.ylim(y_limits[0], y_limits[1])
plt.xlim(-0.1, 0.6)
plt.xlabel('$g_0 - r_0$')
plt.ylabel('$M_g$')
leg = plt.legend(loc='lower center')
for t in leg.get_texts():
    t.set_fontsize('small')
plt.savefig(save_name, papertype='letter')
#plt.show()
plt.close('all')
print '-Done'