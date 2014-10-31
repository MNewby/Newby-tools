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

def do_stuff():
    """ called if main """
    #sky_map()
    #clean_FTO_data()
    #FTO_skyplot(gmin=16.0, gmax=23.5)
    FTO_skyplot(gmin=16.0, gmax=23.5, multistep=0.1)
    #plot_from_file()
    #bhb_tools(out=["sgrplot"])
    #bhb_tomography(step=0.1)
    #giant_tools(out=["c-mag"])
    #lambet_profiles(path="/home/newbym2/Desktop/sdssN-stars/", glb="*bhbs_all*",
    #                lj=0, bj=1, gj=3, lbins=5.0, bbins=1.0, suffix="_bhbs", coords="sgr")
    #lambet_profiles(path="/home/newbym2/Desktop/FTO-stars/", glb="FTO_south_dirty**",
    #                lj=0, bj=1, gj=3, suffix="_TOnew", coords="sgr")


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

def FTO_skyplot(path="/home/newbym2/Desktop/FTO-stars/", gmin=16.0, gmax=22.5, 
                multistep=None):
    """Make a skyplot with updated MSTO stars from SDSS """
    #files = glob.glob(path+"stars*")
    #files.append(path+"South/stars-79-new.txt")
    #files.append(path+"South/stars-82-new.txt")
    #files.append(path+"South/stars-86-new.txt")
    files = glob.glob(path+"FTO_All*")
    files.append(path+"FTO_south_dirty.csv")
    #files = glob.glob(path+"*dirty*")
    out = []
    for f in files:
        #data = np.loadtxt(f, skiprows=1)
        data = np.loadtxt(f, delimiter=",", skiprows=1)
        print "Loaded:", f
        for i in range(data.shape[0]):
            #if data[i,2] < ac.getr(16.0):  continue
            #if data[i,2] > ac.getr(22.5):  continue
            gmag = data[i,3]
            if gmag < gmin:  continue
            if gmag > gmax:  continue
            lam, bet = (ac.lb2sgr(data[i,0], data[i,1], 30.0))[3:5]
            out.append([lam, bet, gmag])
        print "Transformed coordinates:", f
    out = np.array(out)
    if multistep != None:
        runs = np.arange(10.0, 50.0, 5.0)  #for distance
        #runs = np.arange(gmin, gmax, multistep)
        for run in runs:
            print "Starting run", run
            gslice = []
            for i in range(out.shape[0]):
                if ac.getr(out[i,2]) < run:  continue
                if ac.getr(out[i,2]) >= (run+5.0):  continue
                gslice.append(out[i,:])
            gslice = np.array(gslice)
            hist = pp.HistMaker(gslice[:,0], gslice[:,1], 0.5, 0.5, yarea=(-70.0, 40.0), xarea=(20.0, 320.0) )
            hist.savehist(outfile="FTO_"+str(run)+"_hist.csv")
            hist.yflip = 1
            hist.xflip = 1
            hist.ticks = (None, [-60, -40.0, -20.0, 0.0, 20.0, 40.0])
            hist.varea = (0.0, 50.0)
            pp.PlotHist(hist, imfile="FTO_"+str(run)+".png", cbarO='horizontal')
    else:
        hist = pp.HistMaker(out[:,0], out[:,1], 0.5, 0.5, yarea=(-70.0, 40.0), xarea=(20.0, 320.0) )
        hist.savehist(outfile="FTO_All_hist.csv")
        hist.yflip=1
        hist.xflip=1
        hist.labels=(r"$\Lambda$","B")
        hist.ticks = (None, [-40.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0])
        hist.varea=(0.0, 200.0)
        pp.PlotHist(hist, imfile="FTO_All.png", cbarO='horizontal')
    

def plot_from_file(infile="/home/newbym2/Newby-tools/paper-tools/sgr_bhbs.csv"):  #FTO_All_hist.csv
    """ creates a image plot from a file """
    indata = np.loadtxt(infile, delimiter=",")
    #hist = pp.HistMaker([0.0, 1.0], [0.0, 1.0], 0.5, 0.5, yarea=(-70.0, 40.0), xarea=(20.0, 320.0))
    hist = pp.HistMaker([0.0, 1.0], [0.0, 1.0], 1.0, 1.0, yarea=(-70.0, 40.0), xarea=(20.0, 320.0))
    hist.H = indata
    hist.yflip=1
    hist.xflip=1
    hist.labels=(r"$\Lambda$",r"${\rm B}$")
    hist.ticks = (None, [-60, -40, -20, 0, 20, 40])
    #hist.varea=(0.0, 200.0)
    hist.varea=(0.0, 45.0)
    pp.PlotHist(hist, imfile="BHB_All.png", cbarO='horizontal', cax=[0.1, 0.75, 0.8, 0.02])
    
    

def bhb_tools(infile="/home/newbym2/Desktop/sdssN-stars/sgr_bhbs_all.csv", 
                out=[], gmin=15.0, gmax=22.5, suffix="All"):
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
            if data[i,3] > gmax:  continue
            if data[i,3] < gmin:  continue
            sgr.append([data[i,0], data[i,1]])
        sgr = np.array(sgr)
        lam, bet = (ac.lb2sgr(sgr[:,0], sgr[:,1], 30.0))[3:5]
        hist = pp.HistMaker(lam, bet, 1.0, 1.0, xarea=(20.0, 320.0), yarea=(-70.0, 40.0))
        hist.xflip=1
        hist.yflip=1
        hist.varea=(0.0, 6.0)
        hist.savehist("sgr_bhbs.csv")
        #hist.ticks= [[20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0],[-40.0, -20.0, 0.0, 20.0]] 
        # mask the data
        if "mask" in out:
            mask = np.loadtxt("Rhist_sgr.csv", delimiter=",")
            for i in range(len(mask[:,0])):
                for j in range(len(mask[0,:])):
                    if mask[i,j] == 0.0:  hist.H[i,j] = 0.0
        hist.savehist("bhb"+suffix+"_hist.csv")
        pp.PlotHist(hist, imfile="BHB_"+suffix+".png", cbarO='horizontal', cax=[0.1, 0.75, 0.8, 0.02])

def bhb_tomography(minmag=15.0, maxmag=22.5, step=0.5):
    runs = np.arange(minmag, maxmag, step)
    for run in runs:
        print "Starting run", run
        bhb_tools(out=["sgrplot"], gmin=run, gmax=(run+step), suffix=str(run))
    print "### - Done"
    
            
def giant_tools(infile="/home/newbym2/Desktop/sdssN-stars/giant_spec_all.csv", out=[]):
    """ tools for red giants """
    data = np.loadtxt(infile, delimiter=",", skiprows=1)
    ff = 1
    #l,b,dered_g,dered_r,ELODIERVFINAL,ELODIERVFINALERR,loggadop,fehadop,plate,p_ra,p_dec
    if "c-mag" in out:
        gr0 = data[:,2] - data[:,3]
        cmag = pp.HistMaker(gr0,data[:,2], 0.01, 0.1, xarea=(0.2, 1.5), yarea=(16.0, 21.0))
        cmag.yflip=1
        cmag.scale='sqrt'
        cmag.varea=(0.0, 144.0)
        cmag.plot()
    plt.show()
    
def lambet_profiles(path="/home/newbym2/Desktop/sdssN-stars/", glb="*bhbs_all*",
                    lj=0, bj=1, gj=2, lbins=2.5, bbins=0.5, suffix="_bhbs", coords="sgr"):
    """ remixed from plot_profiles in utilities1 """
    files = glob.glob(path+glb)
    print files
    data=[]
    g_cut_low, g_cut_high = 15.0, 22.5
    pb = pr.Progressbar(steps=len(files), prefix="Loading Stars:", suffix=None,
        symbol="#", active="=", brackets="[]", percent=True, size=40)
    for f in files:
        count=0
        stripedata = np.loadtxt(f, delimiter=",", skiprows=1)
        for i in range(stripedata.shape[0]):
            if stripedata[i,gj] < g_cut_low:  continue
            if stripedata[i,gj] > g_cut_high:  continue
            #if ac.SDSS_primary(temp[0],temp[1],wedge,fmt="lb",low=9,high=27)==0:  continue
            data.append([stripedata[i,lj], stripedata[i,bj], stripedata[i,gj] ] )
        pb.updatebar(float(files.index(f)+1)/float(len(files)) )
    data = np.array(data)
    if coords == "sgr":
        count, nStars, pb2 = 0, float(len(data[:,0])), pr.Progressbar(steps=100,
            prefix="Changing Coordinates:", suffix=None, symbol="#", active="=",
            brackets="[]", percent=True, size=40)
        for i in range(data.shape[0]):
            count += 1
            data[i,0], data[i,1] = ac.lb2sgr(data[i,0], data[i,1], 50.0 )[3:5]
            #data[i,2] = data[i,2]
            if count % 1000 == 0:  pb2.updatebar(float(count)/nStars)
    elif coords == "lbr":  pass
        #data[i,0], data[i,1] = ac.lb2GC(data[i,0], data[i,1], 15)
        #data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
    else:  print "!!! Coordinate system not recognized !!!";  sys.exit(2)
    pb2.endbar()
    # loop over Lambda slices and plot Beta histograms
    pb3 = pr.Progressbar(steps=len(files), prefix="Making Slices:", suffix=None,
        symbol="#", active="=", brackets="[]", percent=True, size=40)
    #Ls = (40.0, 130.0, lbins)
    Ls = (200.0, 300.0, lbins)
    Lsteps = int((Ls[1]-Ls[0])/Ls[2])
    for i in range(Lsteps):
        Lname = str(Ls[0]+(Ls[2]*i) )[:6]+suffix
        new = []
        for j in range(len(data[:,0])):
            if data[j,0] < Ls[0]+(Ls[2]*i):  continue
            if data[j,0] > Ls[0]+(Ls[2]*(i+1)):  continue
            new.append([data[j,0], data[j,1], data[j,2]])
        if len(new) < 1:  continue
        new = np.array(new)
        #np.savetxt("stars-lambda-"+("00"+str(i))[-2:]+".txt", new, fmt='%.6f')
        pb3.updatebar(float(i)/float(Lsteps))
        Lhist, Ledges = np.histogram(new[:,1], int(100/bbins), (-60.0, 40.0))
        #output to file
        outstuff = sc.array(zip(Ledges+(bbins/2.0), Lhist))
        np.savetxt(Lname+".out", outstuff, fmt='%.2f')
    pb3.endbar()
    print "Ended Successfully"    

def build_hists(hmin=19.2, hmax=22.5, lambins=5,
    path="/home/newbym2/Dropbox/Research/sgrLetter/FTO_tomography/half_bins/"):
    """ Add together 3-dimensional histograms to make 2-d slices for fitting """
    files = glob.glob(path+"*.csv")
    print files
    for f in files:
        run = f[-8:-4]
        if eval(run) < hmin:  continue
        try: H = H + np.loadtxt(f)
        except NameError:  H = np.loadtxt(f)
    # HERE##########################################################
    nslices = int(H.shape[1]/lambins)
    for i in range(nslices):
        for j in range(lambins):
           pass 
        
            
if __name__ == "__main__":  do_stuff()
