import sys
sys.path.insert(0, '../plotting-tools')
sys.path.insert(0, '../utilities')
import math as ma
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import plot_pack as pp
import astro_coordinates as ac

def make_sim_stream_plot(filename="streamgen_sgr_sim.txt"):
    """ Makes the plot for the simulated streams 
        /home/newbym2/Dropbox/Research/sgrnorth_paper/sgr_separated_stars_MRT.txt"""
    data = np.loadtxt(filename)
    for i in range(len(data[:,0])):
        data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
    sky = pp.HistMaker(data[:,0], data[:,1], xsize=0.5, ysize=0.5, 
        xarea=(100.0, 250.0), yarea=(-10.0, 80.0))
    #sky.varea = (0.0,300.0)
    pp.PlotHist(sky, "sgrsim_radec.png")
    print "Ended Successfully"

if __name__ == "__main__":
    make_sim_stream_plot()
