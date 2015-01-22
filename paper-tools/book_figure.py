import sys
sys.path.insert(0, '../plotting-tools')
sys.path.insert(0, '../utilities')
sys.path.insert(0, '../fitting-tools')
import math as ma
import numpy as np
import scipy as sc
import matplotlib
import matplotlib.pyplot as plt
import plot_pack as pp
import astro_coordinates as ac

def do_stuff():
    #Mgiant_hist()
    book_plot()
    
def Mgiant_hist():
    data = np.loadtxt("/home/newbym/Desktop/Mgiant_wise_sgr.csv", delimiter=",", skiprows=1)
    print data.shape
    #data[:,9] = -1.0*data[:,9]
    #hist = pp.HistMaker(data[:,8], data[:,9], 1.0, 1.0, yarea=(-70.0, 40.0), xarea=(20.0, 320.0) )
    hist = pp.HistMaker(data[:,8], data[:,9], 1.0, 1.0)
    hist.savehist(outfile="Mgiant_hist.csv", fmt='%.1f')
    hist.yflip=1
    hist.xflip=1
    hist.labels=(r"$\Lambda$","B")
    #hist.ticks = (None, [-40.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0])
    hist.varea=(0.0, 200.0)
    pp.PlotHist(hist, imfile="Mgiant_sgr.png", cbarO='horizontal')
        
def book_plot():
    BHB = np.loadtxt("/home/newbym/Desktop/bhbs_all_hist.csv", delimiter=",")
    MSTO = np.loadtxt("/home/newbym/Desktop/FTO_All_hist.csv", delimiter=",")
    RGB = np.loadtxt("/home/newbym/Desktop/Mgiant_hist.csv", delimiter=",")
    r_steps, rx_min, rx_max, ry_min, ry_max = 1.0, 0.0, 360.0, -90.0, 90.0
    f_steps, fx_min, fx_max, fy_min, fy_max = 0.5, 20.0, 320.0, -70.0, 40.0
    b_steps, bx_min, bx_max, by_min, by_max = 1.0, 20.0, 320.0, -70.0, 40.0
    fig = plt.figure(1, figsize=(12,9), dpi=120)
    plt.subplots_adjust(hspace=0.001, wspace=0.001)
    # BHBs
    sp3 = fig.add_subplot(313)
    BHB = BHB[:,::-1]
    #BHB = BHB[::-1,:]
    im3 = sp3.imshow(BHB, cmap=pp.spectral_wb, vmax=45.0)
    xlocs = np.arange(0.0, 301.0, 20.0)
    xlabs = []
    for i in range(len(xlocs)):  xlabs.append("${0:.0f}$".format(xlocs[-1*(i+1)]+20))
    plt.xticks(xlocs, xlabs, fontsize=10)
    ylocs = np.arange(10.0, 120.0, 10.0)
    ylabs = []
    for i in range(len(ylocs)):  
        if i % 2 == 0:  ylabs.append("${0:.0f}$".format(ylocs[i]-70))
        else:  ylabs.append("")
    plt.yticks(ylocs, ylabs, fontsize=10 )
    plt.xlabel(r"$\Lambda$")
    plt.ylabel(r"$B$", rotation=0)
    plt.text(125, 11, "SDSS BHB Selection", fontsize=8, family="serif")
    # RGBs
    sp1 = plt.subplot(311) #, sharex=sp3)
    RGB = RGB[:,::-1]
    RGB = RGB[::-1,:]
    sp1.imshow(RGB[20:131,20:320], cmap=pp.spectral_wb, vmax=30.0)
    plt.xticks(np.arange(20.0, 320.0, 20.0))
    plt.setp(sp1.get_xticklabels(), visible=False)
    ylocs = np.arange(10.0, 120.0, 10.0)
    ylabs = []
    for i in range(len(ylocs)):  
        if i % 2 == 0:  ylabs.append("${0:.0f}$".format(ylocs[i]-70))
        else:  ylabs.append("")
    plt.yticks(ylocs[:-1], ylabs[:-1], fontsize=10 )
    plt.ylabel(r"$B$", rotation=0)
    plt.text(125, 11, "WISE RGB Selection", fontsize=8, family="serif",
        backgroundcolor="w", color="k")
    # MSTOs
    sp2 = plt.subplot(312) #, sharex=sp3)
    MSTO = MSTO[:,::-1]
    #MSTO = MSTO[::-1,:]
    sp2.imshow(MSTO, cmap=pp.spectral_wb, vmax=160.0)
    plt.xticks(np.arange(40.0, 620.0, 40.0))
    plt.setp(sp2.get_xticklabels(), visible=False)
    ylocs = np.arange(20.0, 240.0, 20.0)
    ylabs = []
    for i in range(len(ylocs)):
        if i % 2 == 0:  ylabs.append("${0:.0f}$".format((ylocs[i]/2.0)-70))
        else:  ylabs.append("")
    plt.yticks(ylocs[:-1], ylabs[:-1], fontsize=10 )
    plt.ylabel(r"$B$", rotation=0)
    plt.text(246, 22, "SDSS MSTO Selection", fontsize=8, family="serif")
    plt.savefig("figure_raw.ps")
    plt.savefig("figure_raw.png")
    


if __name__ == "__main__":  do_stuff()
