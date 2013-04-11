import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func
import astro_coordinates as coor
import sdss_visualizers as sv
import poly_fit as poly

'''python script for quickly analyzing data.
Matthew Newby, November 15,2010'''

def do_compare_wedges(file1="stars-82.txt", file2="Stripe82_coadd.csv", stripe=82,
                      mag=0, size=1.0):
    """ Modify if size is not 1.0 """
    one_run = fi.read_data(file1)
    or_l = len(one_run[:,0])
    or_hist = sv.plot_wedge_density(one_run, stripe, q=0.458, r0=19.4,
                                    name="_rho1", mag=mag, plot=0, size=size)
    coadd = fi.read_data(file2)
    ca_l = len(coadd[:,0])
    ca_hist = sv.plot_wedge_density(coadd, stripe, q=0.458, r0=19.4,
                       name="_rho2", mag=mag, plot=0, size=size)
    # Separate into heights
    or_h = or_hist[:,1]
    ca_h = ca_hist[:,1]
    # Divide the first data set by the second
    if len(or_h) < len(ca_h):
        l = len(or_h)
        extra_h = -0.1*sc.ones((len(ca_h)-l))
    else:
        l = len(ca_h)
        extra_h = 0.1*sc.ones((len(or_h)-l))
    diff_h = sc.zeros(l)
    for i in range(l):
        diff_h[i] = ( or_h[i] / ca_h[i] )
    out = sc.zeros((l,3))
    for i in range(l):
        out[i,0], out[i,1] = ca_hist[i,0], diff_h[i]
        out[i,2] = 1.0 #ma.sqrt(or_hist[i,2]*or_hist[i,2] + ca_hist[i,2]*ca_hist[i,2])
    return out

def compare_plot(out, size):
    #Make plot
    fig = plt.figure()
    plt.bar(out[:,0], out[:,1], size)
    #x = sc.arange(0.0, 47.0, 1.0)
    #y = func.SDSS_error_function(x, [1.02, 0.119, -3.53])
    #plt.plot(x, y, 'k-')
    plt.xlabel(r'Sun-centered $r$ (kpc)')
    plt.ylabel(r'$(\rho_{1} / \rho_{2})$ (stars/kpc)')
    #plt.savefig(("diff.ps"), papertype='letter')
    plt.show()
    plt.close('all')

if __name__ == "__main__":
    file1="stars-82.txt"
    file2="stars-82-coadd2.txt"
    file3="matched-82-1arcsec.txt"
    mag = 0
    size = 1.0
    out1 = do_compare_wedges(file1, file3, stripe=82, mag=mag, size=size)
    out2 = do_compare_wedges(file3, file2, stripe=82, mag=mag, size=size)
    out = sc.zeros((len(out2[:,0]),2))
    for i in range(len(out[:,0])):
        out[i,0], out[i,1] = out2[i,0], (out1[i,1]*out2[i,1])
    compare_plot(out, size)
    for i in range(len(out[:,0])):
        out[i,0] = out[i,0] + (size/2.0)
    fi.write_data(out, "eps1_eps2.txt")