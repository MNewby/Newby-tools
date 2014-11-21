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
    #FTO_skyplot(gmin=20.0, gmax=22.5)
    #FTO_skyplot(gmin=16.0, gmax=23.5, multistep=0.1)
    #plot_from_file()
    #bhb_tools(out=["sgrplot"])
    #bhb_tomography(step=0.1)
    #giant_tools(out=["c-mag"])
    #lambet_profiles(path="/home/newbym2/Desktop/sdssN-stars/", glb="*bhbs_all*",
    #                lj=0, bj=1, gj=3, lbins=5.0, bbins=1.0, suffix="_bhbs", coords="sgr")
    #lambet_profiles(path="/home/newbym2/Desktop/FTO-stars/", glb="FTO_south_dirty**",
    #                lj=0, bj=1, gj=3, suffix="_TOnew", coords="sgr")
    #build_hists(hmin=20.0, hmax=22.5, lambins=5, fit_type="double", background="fitline",
    #    path="/home/newbym2/Dropbox/Research/sgrLetter/FTO_tomography/half_bins/",
    #    outpath="/home/newbym2/Dropbox/Research/sgrLetter/plus_15kpc/") #15 kpc
    build_hists(hmin=16.0, hmax=22.5, lambins=5, fit_type="double", background="least",
        lstep=1.0, bstep=1.0,
        path="/home/newbym2/Dropbox/Research/sgrLetter/bhb_tomography/half_bins/",
        outpath="/home/newbym2/Dropbox/Research/sgrLetter/plus_15kpc_bhb/") #15 kpc bhb
    #make_galaxy_plot()

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
        #runs = np.arange(10.0, 50.0, 5.0)  #for distance
        runs = np.arange(gmin, gmax, multistep)
        for run in runs:
            print "Starting run", run
            gslice = []
            for i in range(out.shape[0]):
                if out[i,2] < run:  continue
                if out[i,2] >= (run+multistep):  continue
                gslice.append(out[i,:])
            gslice = np.array(gslice)
            hist = pp.HistMaker(gslice[:,0], gslice[:,1], 0.5, 0.5, yarea=(-70.0, 40.0), xarea=(20.0, 320.0) )
            hist.savehist(outfile="FTO_"+str(run)+"_hist.csv", fmt='%.1f')
            hist.yflip = 1
            hist.xflip = 1
            hist.ticks = (None, [-60, -40.0, -20.0, 0.0, 20.0, 40.0])
            hist.varea = (0.0, 12.0)
            pp.PlotHist(hist, imfile="FTO_"+str(run)+".png", cbarO='horizontal')
    else:
        hist = pp.HistMaker(out[:,0], out[:,1], 0.5, 0.5, yarea=(-70.0, 40.0), xarea=(20.0, 320.0) )
        hist.savehist(outfile="FTO_All_hist.csv", fmt='%.1f')
        hist.yflip=1
        hist.xflip=1
        hist.labels=(r"$\Lambda$","B")
        hist.ticks = (None, [-40.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0])
        hist.varea=(0.0, 200.0)
        pp.PlotHist(hist, imfile="FTO_fits.png", cbarO='horizontal')
    

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
        hist.savehist("sgr_bhbs.csv", fmt='%.1f')
        #hist.ticks= [[20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0],[-40.0, -20.0, 0.0, 20.0]] 
        # mask the data
        if "mask" in out:
            mask = np.loadtxt("Rhist_sgr.csv", delimiter=",")
            for i in range(len(mask[:,0])):
                for j in range(len(mask[0,:])):
                    if mask[i,j] == 0.0:  hist.H[i,j] = 0.0
        hist.savehist("bhb"+suffix+"_hist.csv", fmt='%.1f')
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

def build_hists(hmin=19.2, hmax=22.5, lambins=5, fit_type="double", 
    background="common", lstep=0.5, bstep=0.5,
    path="/home/newbym2/Dropbox/Research/sgrLetter/FTO_tomography/half_bins/",
    outpath="/home/newbym2/Dropbox/Research/sgrLetter/plus_10kpc/"):
    """ Add together 3-dimensional histograms to make 2-d slices for fitting """
    lmin, lmax = 20.0, 320.0
    bmin, bmax = -70.0, 40.0
    files = glob.glob(path+"*.csv")
    print files
    for f in files:
        run = eval(f[-13:-9])
        if run < hmin:  continue
        if run > hmax:  continue
        try: H = H + np.loadtxt(f, delimiter=",")
        except NameError:  H = np.loadtxt(f, delimiter=",")
    nslices = int(H.shape[1]/lambins)
    #pb = pr.Progressbar(steps=H.shape[1], prefix="Compiling Slices:", suffix=None,
    #    symbol="#", active="=", brackets="[]", percent=True, size=40)
    cuts = np.loadtxt(path+"cuts.txt")
    j, hist, lname = 0, np.zeros(H.shape[0], float), "lambda_"+str(50.0)
    for i in range(H.shape[1]):
        lslice = (float(i)*lstep)+lmin
        # Ignore slices with too little data
        if lslice < 50.0:  continue
        if (lslice >= 117.5) and (lslice < 190.0):  continue
        if lslice >= 315.0:  continue
        # Fit slices 
        if j >= lambins:
            np.savetxt(outpath+lname+".txt", hist, fmt='%.2f')
            # FIT STUFF HERE; SKIP IF EMPTY
            if np.sum(hist) > 100:
                xx = np.arange(bmin, bmax+0.00001, bstep)+(bstep*0.5)
                fit_hists(xx, hist, lname, fit_type=fit_type, outpath=outpath, 
                    cuts=cuts, background=background, binsize=bstep)
            lname = "lambda_"+str(lslice)
            hist = np.zeros(H.shape[0], float)
            j = 0
        hist = hist + H[:,i]
        j += 1
        #pb.updatebar(float(i)/float(H.shape[1]))
    #pb.endbar()
    print "### - DONE"


def fit_hists(x_raw, y_raw, name, outfile=None, outpath="", fit_type='double', 
        cuts=None, background="common", binsize=0.5):
    run = eval(name.split("_")[-1])
    # clean data;  don't fit bins with certain criteria
    if cuts != None:
        mask = sc.ones(len(y_raw), dtype=bool)
        for i in range(cuts.shape[0]):
            if cuts[i,0] == run:
                if binsize==0.5:
                    for j in range(cuts[i,1], cuts[i,2]):
                        mask[(j-1)] = cuts[i,3]
                if binsize==1.0:
                    for j in range(cuts[i,1]/2, cuts[i,2]/2):
                        mask[(j-1)] = cuts[i,3]
    else:
        mask = sc.zeros(len(y_raw), dtype=bool)
        for i in range(len(y_raw)):
            if y_raw[i] <= 0.0:  mask[i] = 1  #flag zeros
        #if y_raw[i] > 1000.0:  mask[i] = 1  #flag globular clusters
    # make new array that contains only mask==0
    x, y = [], []
    for i in range(len(y_raw)):
        if mask[i]==0:  x.append(x_raw[i]); y.append(y_raw[i])
    # Set up data
    x, y = np.array(x), np.array(y)
    e = func.poisson_errors(y)
    if background == "common":
        # Make a flat background line set to most common bin height
        hline, hedges = np.histogram(y, bins=40, range=(0.0, 1000.0))
        best = 25.0*np.argmax(hline) - 12.5
        aa, bb = 0.0, best
    elif background == "least":
        ids = np.argsort(y)[:10]
        aa, bb = 0.0, np.mean(y[ids])
    elif background == "fitline":
        # make a line, using the average of 10 bins on each end as the anchor points
        x0, y0, x1, y1 = [], [], [], []
        minbin = 1.0
        if run < 120:  i = (len(y)/2.0) - 5
        else:  i = 1
        while len(x0) < 10:
            if y[i] > minbin:  x0.append(x[i]);  y0.append(y[i])
            i += 1
        if run < 120:  i = (len(y)/2.0) + 5
        elif run > 254:  i = 65
        else:  i = -1
        while len(x1) < 10:
            if y[i] > minbin:  x1.append(x[i]);  y1.append(y[i])
            i -= 1
        xi, yi = np.mean(x0), np.mean(y0)
        xf, yf = np.mean(x1), np.mean(y1)
        #xi, yi = np.mean(x[2:6]), np.mean(y[2:6])
        #xf, yf = np.mean(x[-7:-3]), np.mean(y[-7:-3])
        aa = (yf-yi)/(xf-xi)
        bb = yf - (aa*xf)
    else:
        # null line
        aa, bb = 0.0, 0.0
    #if run < 120:  tf, tm = 200.0, -50.0
    #if run > 120:  tf, tm = 150.0, 40.0
    #  fit it
    fitter = fit.ToFit(x,y,e)
    if fit_type == "double":
        fitter.function=func.double_gauss_line
        fitter.update_params([6.0, 0.0, 5.0, 6.0, -6.0, 5.0, aa, bb])
        fitter.step = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]
        fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma", "slope", "intercept"]
    elif fit_type == "quad":
        fitter.function=func.quad_fat_gauss_line
        fitter.update_params([6.0, 5.0, 1.0, 6.0, -10.0, 1.0, 5.0, 10.0, 5.0, 10.0, aa, bb])
        fitter.step = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
        fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma", "amp","sigma", "amp","sigma","slope","intercept"]
    elif fit_type == "triple":
        fitter.function=func.triple_gauss_floor
        fitter.update_params([6.0, 5.0, 1.0, 6.0, -10.0, 1.0, 2.0, tm, 40.0, tf])
        fitter.step = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 5.0, 10.0]
        fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma", "amp", "mu", "sigma", "floor"]
    print "\n {0}:".format(name)
    print "# - Starting from:", fitter.params, fitter.RR
    print "# - With steps:", fitter.step
    path1 = fit.gradient_descent(fitter, its=10000, line_search=0)
    path2 = fit.MCMC(fitter)
    #errors = fit.get_Hessian_errors(fitter)
    # cut from ToFit
    new_params=fitter.params
    xx = np.arange(-70.0, 40.0, 0.1)
    plt.figure(1)
    plt.bar(x_raw, y_raw, binsize, align='center', color='white', zorder=1)
    plt.bar(fitter.x, fitter.y, binsize, align='center', color='grey', zorder=2)
    plt.plot(xx, fitter.function(xx, new_params), 'k-')
    print new_params[:4]
    plt.plot(xx, func.gaussian_function(xx, new_params[:3]), 'k--')
    plt.plot(xx, func.gaussian_function(xx, new_params[3:6]), 'k--')
    if fit_type == "triple":
        plt.plot(xx, func.gaussian_function(xx, [new_params[6], new_params[7], new_params[8]]), 'k--')
        plt.plot(xx, new_params[9]*sc.ones(len(xx)), "k--")
    if fit_type == "quad":
        plt.plot(xx, func.gaussian_function(xx, [new_params[6], new_params[1], new_params[7]]), 'k--')
        plt.plot(xx, func.gaussian_function(xx, [new_params[8], new_params[4], new_params[9]]), 'k--')
    plt.xlabel("B")
    plt.ylabel("Counts")
    plt.xlim(-70.0, 40.0)
    plt.ylim(0.0, 100.0)
    #plt.text(-65.0, 1100.0, name, fontsize=8 )
    out_name = outpath+name+"_plot"+".png"
    #plt.show()
    plt.savefig(out_name)
    plt.close('all')
    return 1

def make_galaxy_plot():
    """ makes a plot of known globular clusters, dwarf galaxies, and streams"""
    harris_positions = "/home/newbym2/Dropbox/Research/harris_positions.dat"
    hpfile = open(harris_positions, "r")
    gc_names, gc_l, gc_b = [], [], []
    for line in hpfile:
        if line[0]=="#":  continue
        if line.strip()=="":  continue
        gc_names.append(line[1:12].strip())
        gc_l.append(eval(line[52:59]) )
        gc_b.append(eval(line[60:67]) )
    hpfile.close()
    plt.figure()
    for i in range(len(gc_l)):
        if gc_l[i] > 180.0:  gc_l[i] = gc_l[i]-360.0
    plt.scatter(gc_l, gc_b, c="black", s=3, marker="o")
    plt.xlim(-180.0, 180.0)
    plt.ylim(-90.0, 90.0)
    plt.show()
    
    
if __name__ == "__main__":  do_stuff()
