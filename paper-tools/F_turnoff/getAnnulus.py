#! /usr/bin/env python

import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib as mpl
import matplotlib.pyplot as plt
"""This program contains code for building a grey atmospere 'fudge function'
into the Girardi isocrones for ugriz filters.

Matthew Newby (RPI), Aug 20, 2010
"""

def get_annulus(data, c_cut, center, d_cut):
    new_data = []
    for i in range(len(data[:,0])):
        if c_cut != None:
            if ( (data[i,2]-data[i,3]) < c_cut[0] ):  continue
            if ( (data[i,2]-data[i,3]) > c_cut[1] ):  continue
        dist = get_dist((data[i,0], data[i,1]), center)
        if (dist < d_cut[0]):  continue
        if (dist > d_cut[1]): continue
        new_data.append(data[i,:])
    #new_data = sc.array(new_data)
    return len(new_data)
    
def get_dist(point, center):
    # both ra, dec tuples
    x = (point[0] - center[0])*(point[0] - center[0])
    y = (point[1] - center[1])*(point[1] - center[1])
    return np.sqrt(x + y)
    
def make_annulus_series(data, step, limits, center, c_cut=None):
    """ """
    incs = int((limits[1] - limits[0]) / step) + 1
    an_hist = []
    for i in range(incs):
        d_cut = ( (i*step)+limits[0], ((i+1)*step)+limits[0])
        height = get_annulus(data, c_cut, center, d_cut)
        pos = (i*step)+limits[0]
        an_hist.append([pos, height])
    return sc.array(an_hist)

if __name__ == "__main__":
    data = fi.read_data('HR_NGC_5053.csv', ",")
    step = 0.005
    limits = (0.0, 0.14)
    center = (199.109, 17.697)  #Ra, dec
    c_cut = (0.2, 0.5)
    gc_prof = make_annulus_series(data, step, limits, center, c_cut)
    fi.write_data(gc_prof, "annuli_out_med.txt")
    fig = plt.figure(1)
    plt.bar(gc_prof[:,0], gc_prof[:,1], step)
    plt.xlabel("Radius, degrees")
    plt.ylabel("counts")
    plt.show()
    print "# --- Done"