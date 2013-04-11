#! /usr/bin/env python

import scipy as sc
import math as m
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


"""This program Creates a 3D plot of an idealized Milkyway Galaxy
Matthew Newby (RPI), May 11, 2011
"""

#Initial Variables
MW_center = [[7.0], [0.0], [0.0]]
MW_edge = 15.0

def plot_galaxy():
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter([0.0],[0.0],[0.0], c='yellow')
    ax.scatter(MW_center[0], MW_center[1], MW_center[2])
    t = sc.arange(0.0, (2.0*np.pi), 0.01 )
    x = MW_edge * np.cos(t)
    y = MW_edge * np.sin(t)
    z = sc.zeros(len(t))
    ax.plot(x,y,z)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    return 0

def plot_vectors():
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot([0.0,-0.42939],[0.0, 0.89954],[0.0, -0.08031])
    ax.plot([0.0, 0.758242], [0.0, 0.25394], [0.0, -0.6005])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

if __name__ == "__main__":
    plot_galaxy()
    #plot_vectors()