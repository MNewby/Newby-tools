import sys
sys.path.insert(0, '../plotting-tools')
sys.path.insert(0, '../utilities')
sys.path.insert(0, '../fitting-tools')
import math as ma
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import plot_pack as pp
import astro_coordinates as ac
import progress as pr
import glob
import sgr_law as sgr
import fit as fit
import functions as func

def plot_profiles(path="/home/newbym2/Desktop/starfiles"):
    files = glob.glob(path+"/stars*")
    #print files
    data=[]
    pb = pr.Progressbar(steps=len(files), prefix="Loading Stars:", suffix=None,
        symbol="#", active="=", brackets="[]", percent=True, size=40)
    for f in files:
        wedge = int(f.split("-")[1].split(".")[0])
        count=0
        stripedata = open(f, "r")
        for line in stripedata:
            count = count + 1
            if count==1:  continue
            if line.strip()=="":  continue
            temp = line.split()
            for i in range(len(temp)):  temp[i] = float(temp[i])
            if (temp[2] < 30.0) or (temp[2] > 45.0):  continue
            if ac.SDSS_primary(temp[0],temp[1],wedge,fmt="lb",low=9,high=27)==0:  continue
            data.append(temp)
        stripedata.close()
        pb.updatebar(float(files.index(f)+1)/float(len(files)) )
    data = np.array(data)
    count, nStars, pb2 = 0, float(len(data[:,0])), pr.Progressbar(steps=100,
        prefix="Changing Coordinates:", suffix=None, symbol="#", active="=",
        brackets="[]", percent=True, size=40)
    for i in range(len(data[:,0])):
        count = count + 1
        #data[i,0], data[i,1] = ac.lb2GC(data[i,0], data[i,1], 15)
        #data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
        data[i,0], data[i,1] = (ac.lb2sgr(data[i,0], data[i,1], 10.0))[3:5]
        data[i,2] = ac.getg(data[i,2], 4.2)
        if count % 100 == 0:  pb2.updatebar(float(count)/nStars)
    pb2.endbar()
    # loop over Lambda slices and plot Beta histograms
    pb3 = pr.Progressbar(steps=len(files), prefix="Making Slices:", suffix=None,
        symbol="#", active="=", brackets="[]", percent=True, size=40)
    Ls = (200.0, 300.0, 2.5)
    Lsteps = int((Ls[1]-Ls[0])/Ls[2])
    for i in range(Lsteps):
        Lname = str(Ls[0]+(Ls[2]*i) )[:6]
        new = []
        for j in range(len(data[:,0])):
            if data[j,0] < Ls[0]+(Ls[2]*i):  continue
            if data[j,0] > Ls[0]+(Ls[2]*(i+1)):  continue
            new.append(data[j,1])
        new = np.array(new)
        Lhist, Ledges = np.histogram(new, 140, (-30.0, 40.0))
        #output to file
        outstuff = sc.array(zip(Ledges+0.25, Lhist))
        np.savetxt(Lname+".out", outstuff)
        #plot
        plt.figure(1)
        plt.bar(Ledges[:-1], Lhist, 0.5)
        plt.title(Lname)
        plt.xlabel("B")
        plt.savefig(Lname+".png")
        plt.close('all')
        pb3.updatebar(float(i)/float(Lsteps))
    print "Ended Successfully"    
    

def make_sim_stream_plot(filein="stream_50shift.txt", RGB=0, imfile=None):
    """ Makes the plot for the simulated streams
        /home/newbym2/Dropbox/Research/sgrnorth_paper/sgr_separated_stars_MRT.txt"""
    folder = "/home/newbym2/Dropbox/Research/sgrLetter/"
    filename=folder + filein
    #file2="/home/newbym2/Dropbox/Research/sgrnorth_paper/sgr_separated_stars_MRT.txt"
    #file2="streamgen_sgr_sim.txt"
    #file2="streamgen_sgrfidprim.txt"
    file2 = folder+"streamgen_sfp2.txt"
    data = np.loadtxt(filename)
    data2 = []
    datasgr = np.loadtxt(file2)
    for i in range(len(data[:,0])):
        #data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
        #if ac.SDSS_primary(temp[0],temp[1],wedge,fmt="lb",low=9,high=27)==0:  continue
        data[i,0], data[i,1], data[i,2] = (ac.lb2sgr(data[i,0], data[i,1], data[i,2]))[3:]
        lam, bet = new_shift_sgr(data[i,0], data[i,1])  #move 2nd stream to new position
        data2.append([lam, bet, ac.getg(data[i,2], M=4.2)])
    data2 = sc.array(data2)
    for i in range(len(datasgr[:,0])):
        datasgr[i,0], datasgr[i,1] = (ac.lb2sgr(datasgr[i,0], datasgr[i,1], 10.0))[3:5]
        datasgr[i,2] = ac.getg(data[i,2], M=4.2)
    """
    data2 = np.loadtxt(file2)
    for i in range(len(data2[:,0])):
        #data2[i,0], data2[i,1] = ac.lbToEq(data2[i,0], data2[i,1])
        data2[i,0], data2[i,1] = (ac.lb2sgr(data2[i,0], data2[i,1], 10.0))[3:5]
    """
    if RGB==1:  
        data3 = np.concatenate((data2,datasgr), axis=0)
        #data3 = np.concatenate((data2,data), axis=0)
        RGB_plot(data3, imfile=imfile, mask_data="Bhist_sgr.csv", muddle=0)
    else:
        sky = pp.HistMaker(np.concatenate([data[:,0],data2[:,0]]), np.concatenate([data[:,1],data2[:,1]]),
            xsize=0.5, ysize=0.5, xarea=(120.0, 250.0), yarea=(-10.0, 50.0))
        sky.varea = (0.0,200.0)
        sky.H = sky.H + np.random.normal(60.0, 15.0, sky.H.shape)
        #for i in range(sky.H.shape[0]):
        #    if i < 14:   sky.H[i,:] = sky.H[i,:]*0.0; continue
        #    sky.H[i,:] = sky.H[i,:] + 45.0 - (0.5/1.0)*i
        pp.PlotHist(sky, "streamgen_bifExtra_radec.png")
        sky.savehist("streamgen_bifExtra.csv")
    print "Ended Successfully"

def make_total_plot(path="/home/newbym2/Desktop/starfiles", RGB=0, imfile=None):
    files = glob.glob(path+"/stars*")
    #print files
    data=[]
    pb = pr.Progressbar(steps=len(files), prefix="Loading Stars:", suffix=None,
        symbol="#", active="=", brackets="[]", percent=True, size=40)
    for f in files:
        wedge = int(f.split("-")[1].split(".")[0])
        count=0
        stripedata = open(f, "r")
        for line in stripedata:
            count = count + 1
            if count==1:  continue
            if line.strip()=="":  continue
            temp = line.split()
            for i in range(len(temp)):  temp[i] = float(temp[i])
            if ac.SDSS_primary(temp[0],temp[1],wedge,fmt="lb",low=9,high=27)==0:  continue
            data.append(temp)
        stripedata.close()
        pb.updatebar(float(files.index(f)+1)/float(len(files)) )
    data = np.array(data)
    count, nStars, pb2 = 0, float(len(data[:,0])), pr.Progressbar(steps=100,
        prefix="Changing Coordinates:", suffix=None, symbol="#", active="=",
        brackets="[]", percent=True, size=40)
    for i in range(len(data[:,0])):
        count = count + 1
        #data[i,0], data[i,1] = ac.lb2GC(data[i,0], data[i,1], 15)
        #data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
        data[i,0], data[i,1] = (ac.lb2sgr(data[i,0], data[i,1], 10.0))[3:5]
        data[i,2] = ac.getg(data[i,2], 4.2)
        if count % 100 == 0:  pb2.updatebar(float(count)/nStars)
    pb2.endbar()
    if RGB==1:  RGB_plot(data, imfile=imfile)
    else:
        allsky = pp.HistMaker(data[:,0], data[:,1], xsize=0.5, ysize=0.5,
            xarea=(120.0, 250.0), yarea=(-10.0, 50.0))
        #allsky.scale = 'sqrt'
        allsky.varea = (0.0,200.0)
        pp.PlotHist(allsky, "sgrall_GC.png")
        allsky.savehist("SDSSnorthGC.csv")
    print "Ended Successfully"

def make_single_stream_plot(filename):
    data = np.loadtxt(filename)
    for i in range(len(data[:,0])):
        #data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
        #if ac.SDSS_primary(temp[0],temp[1],wedge,fmt="lb",low=9,high=27)==0:  continue
        data[i,0], data[i,1] = (ac.lb2sgr(data[i,0], data[i,1], 10.0))[3:5]
        lam, bet = new_shift_sgr(data[i,0], data[i,1])
        data2.append([lam, bet, data[i,2]])
    # Not done yet!!!
    return -1

def RGB_plot(data, normed=0, imfile=None, mask_data=None, muddle=0):
    xa, ya = (200.0, 300.0), (-40.0, 30.0) #(120.0, 250.0), (-10.0, 50.0)
    #distance cutoffs
    #Wr, Br, Gr, Rr, Kr = 10.0, 20.0, 30.0, 40.0, 50.0
    Wr, Br, Gr, Rr, Kr = 20.0+0.25, 20.66+0.25, 21.33+0.25, 22.0+0.25, 23.0+0.25  #Belokurov 2006 + g-r=0.25
    # bins for each color
    W, B, G, R, K = [], [], [], [], []
    for i in range(len(data[:,0])):
        if   data[i,2] < Wr:  W.append(data[i,:])
        elif data[i,2] < Br:  B.append(data[i,:])
        elif data[i,2] < Gr:  G.append(data[i,:])
        elif data[i,2] < Rr:  R.append(data[i,:])
        else:  pass
    if R == []:  R = [[0.0,0.0,0.0]]  #failsafe
    B = np.array(B)
    Bhist = pp.HistMaker(B[:,0], B[:,1], xsize=0.5, ysize=0.5,
            xarea=xa, yarea=ya)
    if muddle==1:  Bhist.H = Bhist.H + np.absolute(np.random.normal(10.0, 2.0, Bhist.H.shape))
    Bhist.varea = (0.0,200.0)
    Bhist.savehist("Bhist.csv") 
    G = np.array(G)
    Ghist = pp.HistMaker(G[:,0], G[:,1], xsize=0.5, ysize=0.5,
            xarea=xa, yarea=ya)
    if muddle==1:  Ghist.H = Ghist.H + np.absolute(np.random.normal(20.0, 5.0, Ghist.H.shape))
    Ghist.varea = (0.0,200.0)
    Ghist.savehist("Ghist.csv")
    R = np.array(R)
    Rhist = pp.HistMaker(R[:,0], R[:,1], xsize=0.5, ysize=0.5,
            xarea=xa, yarea=ya)
    if muddle==1:  Rhist.H = Rhist.H + np.absolute(np.random.normal(10.0, 5.0, Rhist.H.shape))
    Rhist.varea = (0.0,200.0)
    Rhist.savehist("Rhist.csv")
    if normed == 1:
        Bhist.H = Bhist.H / np.ma.max(Bhist.H)  #Normalize individually
        Ghist.H = Ghist.H / np.ma.max(Ghist.H)
        Rhist.H = Rhist.H / np.ma.max(Rhist.H)
    else:
        norm = max(np.ma.max(Bhist.H), np.ma.max(Ghist.H), np.ma.max(Rhist.H))
        Bhist.H = Bhist.H / norm
        Ghist.H = Ghist.H / norm
        Rhist.H = Rhist.H / norm
    RGB = np.dstack( (Rhist.H, Ghist.H, Bhist.H) )
    # Apply SDSS footprint mask
    if mask_data != None:
        mask = np.loadtxt(mask_data, delimiter=",")  #"SDSSnorth.csv"
        for i in range(len(mask[:,0])):
            for j in range(len(mask[0,:])):
                if mask[i,j] == 0.0:
                    RGB[i,j,:] = 0.0, 0.0, 0.0
    #np.savetxt("RGBout.txt", RGB, delimiter=",")
    plt.figure(1)
    plt.imshow(RGB, origin='lower')
    #xs = np.arange(120, 250, 10)
    #xlocs, xlabels = (xs - 120)*2, []
    #for x in xs:  xlabels.append(str(x))
    #plt.xticks(xlocs, xlabels)
    #ys = np.arange(-10, 51, 10)
    #ylocs, ylabels = (ys + 10)*2, []
    #for y in ys:  ylabels.append(str(y))
    #plt.yticks(ylocs, ylabels)
    if imfile == None:  plt.show()
    else:  plt.savefig(imfile)
    plt.close('all')

def RGB_from_files(mask_data=None, imfile=None):
    suffix = ""
    Rfile = "Rhist"+suffix+".csv"
    Gfile = "Ghist"+suffix+".csv"
    Bfile = "Bhist"+suffix+".csv"
    xa, ya = (200.0, 300.0), (-40.0, 30.0)
    R = np.loadtxt(Rfile, delimiter=",")
    G = np.loadtxt(Gfile, delimiter=",")
    B = np.loadtxt(Bfile, delimiter=",")
    R,G,B = process_image(R,G,B)
    RGB = np.dstack( (R,G,B) )
    # Apply SDSS footprint mask
    if mask_data != None:
        mask = np.loadtxt(mask_data, delimiter=",") #"SDSSnorth.csv"
        for i in range(len(mask[:,0])):
            for j in range(len(mask[0,:])):
                if mask[i,j] == 0.0:
                    RGB[i,j,:] = 0.0, 0.0, 0.0
    #for Sgr;  Sgr coords
    plt.figure(1)
    plt.imshow(RGB) #, origin='lower')
    xs = np.arange(xa[0], xa[1]+1, 20)
    xlocs, xlabels = (xs - xa[0])*2, []
    for x in xs:  xlabels.append(str(x))
    plt.xticks(xlocs, xlabels)
    ys = np.arange(ya[0], ya[1]+1, 20)
    ylocs, ylabels = (ys - ya[0])*2, []
    for y in ys:  ylabels.append(str(y))
    plt.yticks(ylocs, ylabels)
    if imfile == None:  plt.show()
    else:  plt.savefig(imfile)
    plt.close('all')
    
def process_image(R,G,B, stretch=1):
    """processes an RGB image"""
    dialR, dialG, dialB = 1.0, 1.0, 1.0
    gammaR, gammaG, gammaB = 1.0, 1.0, 1.0
    for z in [R,G,B]:
        for i in range(len(z[:,0])):
            for j in range(len(z[0,:])):
                if z[i,j] < 0.2:  z[i,j] = 0.0
    R = (R**gammaR)*dialR
    G = (G**gammaG)*dialG
    B = (B**gammaB)*dialB
    denom = np.ma.max(R) - np.ma.min(R)
    R = (R / denom) - (np.ma.min(R) / denom)
    denom = np.ma.max(G) - np.ma.min(G)
    G = (G / denom) - (np.ma.min(G) / denom)
    denom = np.ma.max(B) - np.ma.min(B)
    B = (B / denom) - (np.ma.min(B) / denom)
    return R, G, B
    


def make_diff_hist():
    """ Subtracts the simulated streams from the real data """
    folder = "/home/newbym2/Dropbox/Research/sgrLetter/"
    data = np.loadtxt("SDSSnorth.csv", delimiter=",")
    sims = ["streamgen_bifGood.csv", "streamgen_bif50Good.csv", "streamgen_bifExtra.csv", "streamgen_bifShift.csv"]
    outs = ["streamgen_ccdiff.png", "streamgen_cc50good.png", "streamgen_ccExtra.png", "streamgen_ccShift.png"]
    #sim = np.loadtxt("simtotbif50m30_radec.csv", delimiter=",")
    #sim = np.loadtxt("simtotbifm2_radec.csv", delimiter=",")
    #sim = np.loadtxt("streamgen_bifExtra.csv", delimiter=",")
    #sim = np.loadtxt("streamgen_bifGood.csv", delimiter=",")
    #sim = np.loadtxt("streamgen_bif50Good.csv", delimiter=",")
    for z in range(len(sims)):
        sim = np.loadtxt(folder+sims[z], delimiter=",")
        for i in range(sim.shape[0]):
            for j in range(sim.shape[1]):
                if data[i,j] < 1.0:  sim[i,j] = 0.0
        sky = pp.HistMaker([1,2], [1,2], xsize=0.5, ysize=0.5,
            xarea=(120.0, 250.0), yarea=(-10.0, 50.0))
        #new = sim #data - sim
        new = data - sim
        sky.H = new
        sky.cmap = 'color'
        sky.varea = (0.0,200.0)
        pp.PlotHist(sky, "sim_"+outs[z])
        #pp.PlotHist(sky, "data.png"); break
    #pp.PlotHist(sky, "streamgen_bif50Good_bwdiff.png")
    print "All done"

def get_bif():
    bifra = [230.0, 220.0, 215.0, 210.0, 200.0, 190.0, 185.0, 180.0, 176.6, 173.3, 170.0, 160.0, 150.0, 140.0, 130.0]
    bifdec = [2.0, 4.0, 6.5, 9.0, 12.0, 15.0, 18.0, 21.0, 23.0, 25.0, 27.0, 28.0, 29.0, 31.0, 31.2]
    stripes = range(9,24)
    for i in range(len(bifra)):
        l,b = ac.EqTolb(bifra[i], bifdec[i])
        for s in stripes:
            if test_primary(l,b,s,low=9,high=23) == 1:
                mu, nu = ac.lb2GC(l,b,s)
                print "RA {0}, dec {1}, is in stripe {2}; mu {3}, nu {4}".format(bifra[i], bifdec[i], s, mu, nu)


def get_sgr_curves():
    #input in pixels
    primary = np.array([[5.0, 4.4], [9.7, 3.6], [14.4, 1.35]])
    secondary = np.array([[3.3,6.5],[10.4,5.0],[14.4,2.0]])
    #convert pixels to ra/dec
    for data in [primary, secondary]:
        data[:,0] = data[:,0]*0.5*12.5 + 120.0
        data[:,1] = data[:,1]*0.5*12.5 - 10.0
        a, b = np.zeros((3,3)),  data[:,1]
        for i in range(3):
            a[i,:] = data[i,0]*data[i,0], data[i,0], 1.0
        print a
        print np.linalg.solve(a,b)
    print primary
    print secondary
    # decl as a function of RA:  (decl = a*RA*RA + b*RA + c)
    prime = [ -5.25124491e-03, 1.57254414e+00, -1.00216868e+02]
    second = [ -7.76551199e-03, 2.31737724e+00, -1.41690141e+02]
    t = np.arange(120.0, 250.0, 1.0)
    u1 = t*t*prime[0] + t*prime[1] + prime[2]
    u2 = t*t*second[0] + t*second[1] + second[2]
    plt.figure()
    plt.plot(t,u1)
    plt.plot(t, u2)
    plt.ylim(-5.0, 50.0)
    plt.show()
    plt.close('all')

def shift_sgr(filein="streamgen_bif50Good.txt", fileout="stream_50shift.txt"):
    prime = [ -5.25124491e-03, 1.57254414e+00, -1.00216868e+02]
    second = [ -7.76551199e-03, 2.31737724e+00, -1.41690141e+02]
    #file1="/home/newbym2/Dropbox/Research/sgrLetter/"+filein
    file1="/home/newbym2/Dropbox/Research/sgrLetter/"+filein
    data = np.loadtxt(file1)
    for i in range(len(data[:,0])):
        data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
        deldel = data[i,1] - (data[i,0]*data[i,0]*prime[0] + data[i,0]*prime[1] + prime[2])
        data[i,1] = deldel + (data[i,0]*data[i,0]*second[0] + data[i,0]*second[1] + second[2])
        data[i,0], data[i,1] = ac.EqTolb(data[i,0], data[i,1])
    np.savetxt(fileout, data, delimiter=" ")

def new_shift_sgr(lam, bet):
    """ shifts Sgr according in Law/Majewski coords """
    cm = (0.2*lam) - 40.0
    cm = (-0.12*cm) + 3.2
    lam = lam + (cm*5.0)
    bet = bet - 10.0
    return lam, bet
    

def batch_shift():
    infiles = ["streamgen_sgr_sim.txt", "streamgen_sgrwide.txt", "streamgen_sgrsmall.txt", "streamgen_sgrbig.txt"]
    outfiles = ["stream_shift.txt", "stream_shiftwide.txt", "stream_shiftsmall.txt", "stream_shiftbig.txt"]
    for i in range(len(infiles)):
        infile, outfile = infiles[i], outfiles[i]
        shift_sgr(infile, outfile)
        make_sim_stream_plot(filein=outfile, RGB=1, imfile=outfile[:-4]+".png")
        print "Done with {0}, {1}".format(infile, outfile)


def sgr_rv():
    data = np.loadtxt("/home/newbym2/Dropbox/Research/sgrLetter/sgr_spec.csv", delimiter=",")
    #dered_g,dered_r,l,b,ELODIERVFINAL,ELODIERVFINALERR
    g0, r0 = data[:,0], data[:,1]
    l, b = data[:,2], data[:,3]
    rv, rv_err = data[:,4], data[:,5]
    d = ac.getr(g0)
    # Transform to vgsr from Yanny+ 2009
    vgsr = ac.rv_to_vgsr(rv,l,b)
    X,Y,Z, lsgr, bsgr, r_sgr = ac.lb2sgr(l, b, d)
    for w in [0.5, 1.0, 2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0, 30.0, 50.0]:
        keep = []
        for i in range(len(data[:,0])):
            if abs(bsgr[i]) < w:  
                if g0[i] < 20.0:  continue
                if g0[i] > 23.0:  continue
                keep.append(vgsr[i])
        hist, edges = np.histogram(np.array(keep), bins=60, range=(-300.0, 300.0))
        y, x = hist, edges[:-1] 
        e = func.poisson_errors(y)
        fitter = fit.ToFit(x,y,e)
        fitter.function=func.double_gaussian_one_fixed
        fitter.update_params([3.0, -120.0, 30.0, 2.0, 0.0, 120.0])
        fitter.step = [1.0, 10.0, 1.0, 1.0, 0.0, 0.0]
        fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma"]
        path1 = fit.gradient_descent(fitter, its=10000, line_search=0)
        path2 = fit.MCMC(fitter)
        new_params=fitter.params
        xx = sc.arange(-300.0, 300.0, 1.0)
        yy = func.double_gaussian_one_fixed(xx, new_params)
        y1 = func.gaussian_function(xx, new_params[:3])
        y2 = func.gaussian_function(xx, new_params[3:])
        fig = plt.figure()
        plt.bar(edges[:-1], hist, width=10.0, color="white")
        plt.plot(xx,yy, "k-")
        plt.plot(xx,y1, "k--")
        plt.plot(xx,y2, "k--")
        plt.title("cut width = "+str(w))
        plt.xlabel(r"$v_{\rm gsr}$", fontsize=16)
        plt.ylabel(r"Counts")
        #plt.ylim(0.0, 60.0)
        #plt.show()
        plt.savefig("/home/newbym2/Dropbox/Research/sgrLetter/sgr_spec/r_cut_relative"+str(w)+".png")
        plt.close()    
    
    
def proj_test():
    #data = np.loadtxt("/home/newbym2/Dropbox/Research/sgrLetter/sgr_spec.csv", delimiter=",")
    x = sc.arange(0.0, 359.0, 0.1)
    y = sc.zeros(len(x))
    y2 = sc.arange(-90.0, 90.0, 0.1)
    x2 = sc.zeros(len(y2))
    y3 = sc.arange(-90.0, 90.0, 0.1)
    x3 = 180.0*sc.ones(len(y3))
    y4 = sc.arange(-90.0, 90.0, 0.1)
    x4 = 359.0*sc.ones(len(y4))
    newx = np.append(x,x2)
    newx = np.append(newx, x3)
    newx = np.append(newx, x4)
    newy = np.append(y,y2)
    newy = np.append(newy, y3)
    newy = np.append(newy, y4)
    fig = plt.figure()
    #xx, yy = ac.equirec_projection(newx-180.0, newy)
    xx, yy = ac.equal_area_projection(newx-180.0, newy)
    #xx, yy = ac.tripel_projection(newx-180.0, newy)
    #xx, yy = ac.sin_projection(newx-180.0, newy)
    plt.scatter(xx, yy, marker="o", s=1, c="k")
    plt.show()
    plt.close()
    

if __name__ == "__main__":
    #shift_sgr(filein="streamgen_sgrfidprim.txt", fileout="stream_shiftfid.txt")
    #make_sim_stream_plot(filein="streamgen_sfp_bigish.txt", RGB=1) #, imfile="sgr_new.png")
    #make_total_plot(RGB=1)
    #make_diff_hist()
    #get_bif()
    #get_sgr_curves()
    #batch_shift()
    #RGB_from_files(mask_data="Rhist_sgr.csv", imfile="new.png")
    #plot_profiles()    
    #proj_test()
    sgr_rv()
    
"""
RA 230.0, dec 2.0, is in stripe 11; mu 229.965574947, nu 0.23156178591
RA 220.0, dec 4.0, is in stripe 12; mu 219.902390269, nu -0.0990566496675
RA 215.0, dec 6.5, is in stripe 13; mu 214.787585638, nu -0.00444745570075
RA 210.0, dec 9.0, is in stripe 14; mu 209.671526225, nu -0.0792531928135
RA 200.0, dec 12.0, is in stripe 15; mu 199.664792156, nu -0.0866711273422
RA 190.0, dec 15.0, is in stripe 16; mu 189.829219533, nu 0.0545069430007
RA 185.0, dec 18.0, is in stripe 17; mu 185.0, nu 0.5
RA 180.0, dec 21.0, is in stripe 18; mu 180.332045544, nu 1.06962833929
RA 176.6, dec 23.0, is in stripe 19; mu 177.271410366, nu 0.716530158374
RA 173.3, dec 25.0, is in stripe 20; mu 174.409199595, nu 0.455971801713
RA 170.0, dec 27.0, is in stripe 21; mu 171.666854115, nu 0.303227642608
RA 160.0, dec 28.0, is in stripe 22; mu 163.089514843, nu 0.370314382846
RA 150.0, dec 29.0, is in stripe 23; mu 154.880506976, nu 1.37166887545
RA 140.0, dec 31.0, is in stripe 23; mu 147.430539542, nu 6.24134308136
RA 130.0, dec 31.2, is in stripe 23; mu 139.648046388, nu 9.97924011117
"""
