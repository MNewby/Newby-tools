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
import progress as pr
import glob
import sgr_law as sgr
import fit as fit
import functions as func


def sky_map(infile="/home/newbym2/Desktop/sdssN-stars/sdssN_sky.csv"):
    """Map a variable in Lambda, Beta  """
    data = np.loadtxt(infile, delimiter=",", skiprows=1)
    l, b = data[:,0], data[:,1]
    g_ext, r_ext = data[:,4], data[:,5]
    g_sky, r_sky = data[:,6], data[:,7]
    print np.ma.min(g_ext), np.ma.max(g_ext)
    lam, bet = ac.lb2sgr(l, b, 30.0)[3:5]
    plt.figure()
    sc = plt.scatter(lam, bet, c=r_ext, cmap=pp.spectral_wb, edgecolors='none', 
                    s=3, alpha=0.3, vmin=0, vmax=0.4)
    plt.colorbar(sc)
    plt.xlim(230.0, 255.0); plt.ylim(15.0, -10.0)
    plt.show()

def clean_FTO_data(path="/home/newbym2/Desktop/FTO-stars/"):
    """ compactify FTO data """
    files = glob.glob(path+"FTO_All*")
    print files
    out = []
    for f in files:
        data = np.loadtxt(f, delimiter=",", skiprows=1)
        for i in range(data.shape[0]):
            out.append([data[i,0], data[i,1], data[i,3] ])
        print "Completed", f
    out = np.array(out)
    np.savetxt(path+"stars-MSTO-N-all.txt", out, fmt='%.5f')
    print "### - DONE"

def bhb_tools(infile="/home/newbym2/Desktop/sdssN-stars/sgr_bhbs_all.csv", out=[]):
    """ bhb stuff """
    data = np.loadtxt(infile, delimiter=",", skiprows=1)
    ff = 1
    if "c-c" in out:
        ug0 = data[:,2] - data[:,3]
        gr0 = data[:,3] - data[:,4]
        plt.figure(ff);  ff+=1
        plt.scatter(gr0, ug0, c="k", s=2)
        plt.xlabel(r"$(g-r)_0$");  plt.ylabel(r"$(u-g)_0$")
        plt.xlim(-0.5, 1.5);  plt.ylim(3.0, 0.0)
    if "c-mag" in out:
        #rr = ac.getr(data[:,3], 0.7)  ## BHB M  #Not working right?
        gr0 = data[:,3] - data[:,4]
        cmag = pp.HistMaker(gr0,data[:,3], 0.01, 0.1)
        cmag.yflip=1
        cmag.scale='sqrt'
        cmag.varea=(0.0, 2000.0)
        cmag.plot()
    plt.show()
    if "sgrplot" in out:
        # use only primaries
        sgr = []
        for i in range(data.shape[0]):
            #if ac.SDSS_primary(data[i,0],data[i,1],wedge,fmt="lb",low=9,high=27)==0:  continue
            if data[i,3] > 22.5:  continue
            if data[i,3] < 15.0:  continue
            sgr.append([data[i,0], data[i,1]])
        sgr = np.array(sgr)
        lam, bet = (ac.lb2sgr(sgr[:,0], sgr[:,1], 30.0))[3:5]
        hist = pp.HistMaker(lam, bet, 1.0, 1.0, xarea=(20.0, 320.0), yarea=(-70.0, 40.0))
        hist.yflip=1
        hist.varea=(0.0, 50.0)
        #hist.ticks= [[20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0],[-40.0, -20.0, 0.0, 20.0]] 
        # mask the data
        if "mask" in out:
            mask = np.loadtxt("Rhist_sgr.csv", delimiter=",")
            for i in range(len(mask[:,0])):
                for j in range(len(mask[0,:])):
                    if mask[i,j] == 0.0:  hist.H[i,j] = 0.0
        hist.plot()
            


def sgrN_dr8(path="/home/newbym2/Desktop/FTO-stars"):
    """ """
    pass


            
if __name__ == "__main__":
    #sky_map()
    #clean_FTO_data()
    bhb_tools(out=["sgrplot"])
