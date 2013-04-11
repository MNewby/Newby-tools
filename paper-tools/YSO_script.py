import math as ma
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import sys
import files as fi

"""Script containing code for the YSO-SDSS project
    Matthew Newby, Dec 30, 2011"""

'''
l, b,
psfMag_u, psfMag_g,psfMag_r,psfMag_i,psfMag_z,
extinction_u,extinction_g,extinction_r,extinction_i,extinction_z
'''

ug_limits, gr_limits = (-0.5, 3.5, 0.02), (-0.5, 2.0, 0.02)
ug_bins = int((ug_limits[1]-ug_limits[0])/ug_limits[2]) + 1
gr_bins = int((gr_limits[1]-gr_limits[0])/gr_limits[2]) + 1
    
def make_hists(filename):
    # 360 degree histogram in l
    l_hist = sc.zeros(360, float)
    # u-g from -0.5 to 3.5, in 0.1 increments; g-r from -0.5 to 2.0, in 0.1 increments
    dered_hist = sc.zeros((ug_bins, gr_bins), float)
    # same hist, but for un-dereddened stars
    red_hist = sc.zeros((ug_bins, gr_bins), float)
    # Get data
    counter = 0
    readfile = open(filename, "r")
    for line in readfile:
        counter = counter + 1
        if counter % 100000 == 0:  print "On line {0}".format(counter)
        if line.strip() == "":  continue
        if line.strip()[0] == "#":  continue
        holder = line.split(",")
        """ Get l """
        l = eval(holder[0])
        while l > 360.0:
            l = l - 360.0
        while l < 0.0:
            l = l + 360.0
        lindex = int(l)
        l_hist[lindex] = l_hist[lindex] + 1.0
        """ Fill color-color plots """
        ug = eval(holder[2])-eval(holder[3])
        gr = eval(holder[3])-eval(holder[4])
        ug_red = eval(holder[7])-eval(holder[8])
        gr_red = eval(holder[8])-eval(holder[9])
        ug_dered = ug - ug_red
        gr_dered = gr - gr_red
        # undereddened plot
        ugi = int((ug - ug_limits[0])/ug_limits[2])
        if (ugi < 0) or (ugi > ug_bins-1):  continue
        gri = int((gr - gr_limits[0])/gr_limits[2])
        if (gri < 0) or (gri > gr_bins-1):  continue
        red_hist[ugi, gri] = red_hist[ugi, gri] + 1.0 #if this fails, reverse the axes
        # dereddened plot
        ugj = int((ug_dered - ug_limits[0])/ug_limits[2])
        if (ugj < 0) or (ugj > ug_bins-1):  continue
        grj = int((gr_dered - gr_limits[0])/gr_limits[2])
        if (grj < 0) or (grj > gr_bins-1):  continue
        dered_hist[ugj, grj] = dered_hist[ugj, grj] + 1.0 #if this fails, reverse the axes
    readfile.close()
    return l_hist, red_hist, dered_hist

def save_hists(l_hist=None, red_hist=None, dered_hist=None):
    if l_hist != None:
        print "# --- length of l_hist: {0}".format(len(l_hist))
        l_out = sc.zeros((360, 2),float)
        for i in range(len(l_out[:,0])):
            l_out[i,0] = float(i)
            l_out[i,1] = l_hist[i]
        fi.write_data(l_out, "l_hist_data.csv", delimiter=",",
                      header="angle, number of stars")
    if red_hist != None:
        print "# --- shape of red_hist: {0}".format(red_hist.shape)
        fi.write_data(red_hist, "red_hist_data.csv", delimiter=",")
    if dered_hist != None:
        print "# --- shape of dered_hist: {0}".format(dered_hist.shape)
        fi.write_data(dered_hist, "dered_hist_data.csv", delimiter=",")

def make_lhist_plot(l_hist):
    fig = plt.figure(1)
    if len(l_hist[0,:]) < 2:
        x = sc.arange(0.0, 361.0, 1.0)
        y = l_hist
    else:
        x, y = l_hist[:,0], l_hist[:,1]
    plt.bar(x, y, width=1.0, ec='k', fc=None, fill=False) #, hatch='/')
    plt.xlabel("Galactic longitude, l")
    print plt.xticks(), len(x), len(y)
    plt.xlim(0.0, 360.0)
    plt.savefig("YSO_l_hist.ps",papertype='letter')
    plt.close('all')

#ug_limits, gr_limits = (-0.5, 3.5, 0.02), (-0.5, 2.0, 0.02)
#ugi = int((ug - ug_limits[0])/ug_limits[2])
#gri = int((gr - gr_limits[0])/gr_limits[2])

def make_cc_plots(dered_hist, red_hist):
    plt.figure(1)
    plt.imshow(dered_hist, interpolation='nearest') #, cmap='gray')
    """ Set up g-r (x) axis """
    plt.xlabel(r"$(g-r)_0$")
    xplace = sc.arange(0.0, gr_bins, 20.0)
    xread = []
    for i in range(len(xplace)):
        xread.append(str(round( ((xplace[i]*gr_limits[2]) + gr_limits[0]), 1)))
    plt.xticks(xplace, xread)
    """ Set up u-g (y) axis """
    plt.ylabel(r"$(u-g)_0$")
    yplace = sc.arange(0.0, ug_bins, 20.0)
    yread = []
    for i in range(len(yplace)):
        yread.append(str(round( ((yplace[i]*ug_limits[2]) + ug_limits[0]), 1)))
    plt.yticks(yplace, yread)
    #plt.savefig("YSO_dered_hist.ps",papertype='letter')
    """ ---- SECOND PLOT  ---- """
    plt.figure(2)
    plt.imshow(red_hist, interpolation='nearest') #, cmap='gray')
    plt.xlabel(r"$g-r$")
    plt.xticks(xplace, xread)
    plt.ylabel(r"$u-g$")
    plt.yticks(yplace, yread)
    #plt.savefig("YSO_red_hist.ps",papertype='letter')
    plt.show()
    plt.close('all')

def unload_hists(l_in=None, red_in=None, dered_in=None):
    if l_in != None:  l_hist = fi.read_data(l_in, ",")
    else:  l_hist = None
    if red_in != None:  red_hist = fi.read_data(red_in, ",")
    else:  red_hist = None
    if dered_in != None:  dered_hist = fi.read_data(dered_in, ",")
    else:  dered_hist = None
    return l_hist, red_hist, dered_hist

if __name__ == "__main__":
    #l, red, dered = make_hists(sys.argv[1])
    #save_hists(l, red, dered)
    l, red, dered = unload_hists("l_hist_data.csv", "red_hist_data.csv", "dered_hist_data.csv")
    #make_lhist_plot(l)
    make_cc_plots(dered, red)
