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

def plot_profiles(path="/home/newbym2/Desktop/starfiles", savewedge=False):
    files = glob.glob(path+"/stars*")
    #print files
    data=[]
    r_cut_low, r_cut_high = ac.getr(16.0), ac.getr(22.5) #30.0, 45.0
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
            if (temp[2] < r_cut_low) or (temp[2] > r_cut_high):  continue
            if ac.SDSS_primary(temp[0],temp[1],wedge,fmt="lb",low=9,high=27)==0:  continue
            data.append(temp)
        stripedata.close()
        pb.updatebar(float(files.index(f)+1)/float(len(files)) )
    datar = np.array(data)
    data = np.zeros(datar.shape, float)
    count, nStars, pb2 = 0, float(len(data[:,0])), pr.Progressbar(steps=100,
        prefix="Changing Coordinates:", suffix=None, symbol="#", active="=",
        brackets="[]", percent=True, size=40)
    for i in range(len(data[:,0])):
        count = count + 1
        #data[i,0], data[i,1] = ac.lb2GC(data[i,0], data[i,1], 15)
        #data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
        data[i,0], data[i,1] = (ac.lb2sgr(datar[i,0], datar[i,1], 10.0))[3:5]
        data[i,2] = ac.getg(datar[i,2], 4.2)
        if count % 100 == 0:  pb2.updatebar(float(count)/nStars)
    pb2.endbar()
    # loop over Lambda slices and plot Beta histograms
    pb3 = pr.Progressbar(steps=len(files), prefix="Making Slices:", suffix=None,
        symbol="#", active="=", brackets="[]", percent=True, size=40)
    Ls = (200.0, 300.0, 2.5)
    Lsteps = int((Ls[1]-Ls[0])/Ls[2])
    for i in range(Lsteps):
        Lname = str(Ls[0]+(Ls[2]*i) )[:6]+"_NEW"
        new = []
        for j in range(len(data[:,0])):
            if data[j,0] < Ls[0]+(Ls[2]*i):  continue
            if data[j,0] > Ls[0]+(Ls[2]*(i+1)):  continue
            if savewedge:  new.append([datar[j,0], datar[j,1], data[j,2]])
            else:  new.append(data[j,1])
        new = np.array(new)
        if savewedge:
            np.savetxt("stars-lambda-"+("00"+str(i))[-2:]+".txt", new, fmt='%.6f')
        else:
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

def make_total_plot(path="/home/newbym2/Desktop/starfiles", RGB=0, imfile=None, 
                    rcut=None, outfile=None):
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
            if rcut != None:
                if temp[2] < rcut[0]:  continue
                if temp[2] > rcut[1]:  continue
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
            xarea=(200.0, 300.0), yarea=(-40.0, 30.0))
            #xarea=(120.0, 250.0), yarea=(-10.0, 50.0))
        #allsky.scale = 'sqrt'
        allsky.cmap="bw"
        allsky.yflip = 1
        allsky.varea = (0.0,100.0)
        if outfile==None:  pp.PlotHist(allsky, "sgrall_GC.png")
        else:  pp.PlotHist(allsky, outfile)
        #allsky.savehist("SDSSnorthGC.csv")
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
    #Wr, Br, Gr, Rr, Kr = 20.0+0.25, 20.66+0.25, 21.33+0.25, 22.0+0.25, 23.0+0.25  #Belokurov 2006 + g-r=0.25
    Wr, Br, Gr, Rr, Kr = 20.0, 20.1, 20.4, 20.5, 50.0 #spec cut limits
    # bins for each color
    W, B, G, R, K = [], [], [], [], []
    for i in range(len(data[:,0])):
        if   data[i,2] < Wr:  W.append(data[i,:])
        elif data[i,2] < Br:  B.append(data[i,:])
        elif data[i,2] < Gr:  G.append(data[i,:])
        elif data[i,2] < Rr:  R.append(data[i,:])
        else:  pass
    if R == []:  R = [[0.0,0.0,0.0]]  #failsafe
    if G == []:  G = [[0.0,0.0,0.0]]  #failsafe
    if B == []:  B = [[0.0,0.0,0.0]]  #failsafe
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

def RGB_from_files(mask_data=None, imfile=None, fitfile=False):
    wd = "/home/newbym2/Dropbox/Research/sgrLetter/fit_results/"
    fname = "smarter_slide_back_quad.txt"
    suffix = "_sgr2"
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
    plt.figure(1, figsize=(12,9), dpi=100)
    plt.imshow(RGB, zorder=1) #, origin='lower')
    # lambda slices to skip primary or secondary stream plots
    #skip1 = [295.0]  #double
    #skip2 = [210.0]  #double
    skip1 = [297.5, 295.0, 292.5, 290.0] #quad:  
    skip2 = [200.0, 202.5, 212.5, 297.5] #quad:
    if fitfile:
        fitdata = np.loadtxt(wd+fname, delimiter=",")
        lams = (fitdata[:,0] - 200.0) / 0.5
        bet1 = (fitdata[:,2] + 40.0) / 0.5
        bet1e = fitdata[:,3] / 0.5
        bet1ee = fitdata[:,8] / 0.5
        bet2 = (fitdata[:,5] + 40.0) / 0.5
        bet2e = fitdata[:,6] / 0.5
        bet2ee = fitdata[:,10] / 0.5
        # Elminate fits with small amplitudes
        for i in range(len(lams)):
            if (abs(fitdata[i,1]) > 5.0) and (fitdata[i,0] not in skip1):  
                plt.errorbar(lams[i], bet1[i], yerr=bet1e[i], ecolor="white", marker="o", 
                     mec='black', mfc='white', ms=2, ls=" ", mew=0.5, zorder=3)
            if (abs(fitdata[i,4]) > 5.0) and (fitdata[i,0] not in skip2):
                plt.errorbar(lams[i], bet2[i], yerr=bet2e[i], ecolor="white", marker="o", 
                     mec='black', mfc='white', ms=2, ls=" ", mew=0.5, zorder=3)
            if fname[-8:-4] == "quad":
                if (abs(fitdata[i,7]) > 5.0) and (fitdata[i,0] not in skip1):
                    plt.errorbar(lams[i], bet1[i], yerr=bet1ee[i], ecolor="red", fmt=None, capsize=10.0, zorder=2)
                if (abs(fitdata[i,9]) > 5.0) and (fitdata[i,0] not in skip2):
                    plt.errorbar(lams[i], bet2[i], yerr=bet2ee[i], ecolor="magenta", fmt=None, capsize=10.0, zorder=2) 
    xs = np.arange(xa[0], xa[1]+1, 10)
    xlocs, xlabels = (xs - xa[0])*2, []
    for x in xs:  xlabels.append(str(int(x)))
    plt.xticks(xlocs, xlabels, fontsize=15)
    ys = np.arange(ya[0], ya[1]+1, 10)
    ylocs, ylabels = (ys - ya[0])*2, []
    for y in ys:  ylabels.append(str(int(y)))
    plt.yticks(ylocs, ylabels, fontsize=15)
    plt.xlabel(r"$\Lambda$", fontsize=20)
    plt.ylabel(r"$B$", rotation='horizontal', fontsize=20)
    plt.xlim(0.0, 200.0)
    plt.ylim(140.0, 0.0)
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


def sgr_rv_fits():
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
    
def sgr_rv_skyplot():
    data = np.loadtxt("/home/newbym2/Dropbox/Research/sgrLetter/sgr_spec_all.csv", delimiter=",")
    #dered_g,dered_r,l,b,ELODIERVFINAL,ELODIERVFINALERR,plate,p_ra,p_dec
    g0, r0 = data[:,0], data[:,1]
    l, b = data[:,2], data[:,3]
    rv, rv_err = data[:,4], data[:,5]
    plates = data[:,6]
    p_ra, p_dec = data[:,7], data[:,8]
    d = ac.getr(g0)
    # Transform to vgsr from Yanny+ 2009
    vgsr = ac.rv_to_vgsr(rv,l,b)
    X,Y,Z, lsgr, bsgr, r_sgr = ac.lb2sgr(l, b, d)
    pl,pb = ac.EqTolb(p_ra, p_dec)
    pd = 20.0*sc.ones(len(pl))
    pX,pY,pZ, plsgr, pbsgr, pr_sgr = ac.lb2sgr(pl, pb, pd)
    mean, stdev, nsig = -118.233333, 30.222222, 1.0
    #plates = plate_center_utils()
    locs = []
    for i in range(len(data[:,0])):
        if g0[i] < 20.0:  continue
        if g0[i] > 23.0:  continue
        #if lsgr[i] > 230.0:  continue
        locs.append([lsgr[i], bsgr[i], vgsr[i], plsgr[i], pbsgr[i] ])
    locs = np.array(locs)
    sgr = []
    for i in range(len(locs[:,0])):
        if locs[i,2] < (mean - stdev*nsig):  continue
        if locs[i,2] > (mean + stdev*nsig):  continue
        sgr.append([locs[i,0], locs[i,1]])
    sgr = np.array(sgr)
    plt.figure()
    plt.scatter(locs[:,0], locs[:,1], c="k", s=1)
    plt.scatter(sgr[:,0], sgr[:,1], c='r', s=10)
    plt.scatter(locs[:,3], locs[:,4], c='b', s=10)
    plt.plot([190.0, 250.0, 250.0, 190.0],[10.0, 10.0, -10.0, -10.0], "b-")
    plt.plot([190.0, 320.0], [0.0, 0.0], "b--")
    plt.text(195.0, 33.0, "$\mu=$-118.2, $\sigma=$30.2\n{0} stars in {1}$\sigma$".format(len(sgr), nsig), fontsize=16)
    plt.xlim(190.0, 320.0)
    #plt.ylim(-75.0, 40.0)
    plt.ylim(40.0, -75.0)
    plt.xlabel(r"$\Lambda$", fontsize=16)
    plt.ylabel(r"$B$", fontsize=16)
    plt.show()
    """
    hist, edges = np.histogram(sgr[:,1], bins=23, range=(-75.0, 40.0))
    plt.bar(edges[:-1], hist, width=5.0, color="grey")
    plt.xlabel("B", fontsize=16)
    plt.text(-50.0, 25.0, "Stars within "+str(nsig)+r"$\sigma$ of mean $v_{\rm gsr}$; $\Lambda<230$", fontsize=12)
    plt.show()"""
    plt.close()
    

def sgr_rv_cut():
    data = np.loadtxt("/home/newbym2/Dropbox/Research/sgrLetter/sgr_spec_all.csv", delimiter=",")
    #dered_g,dered_r,l,b,ELODIERVFINAL,ELODIERVFINALERR,plate,p_ra,p_dec
    g0, r0 = data[:,0], data[:,1]
    l, b = data[:,2], data[:,3]
    rv, rv_err = data[:,4], data[:,5]
    plates = data[:,6]
    p_ra, p_dec = data[:,7], data[:,8]
    d = ac.getr(g0)
    # Transform to vgsr from Yanny+ 2009
    vgsr = ac.rv_to_vgsr(rv,l,b)
    X,Y,Z, lsgr, bsgr, r_sgr = ac.lb2sgr(l, b, d)
    mean, stdev, nsig = -118.233333, 30.222222, 2.0
    keep = []
    for i in range(len(data[:,0])):
        if abs(bsgr[i]) < 10.0:  
            if g0[i] < 20.0:  continue
            if g0[i] > 23.0:  continue
            if vgsr[i] < (mean - stdev*nsig):  continue
            if vgsr[i] > (mean + stdev*nsig):  continue
            keep.append([(g0[i]-r0[i]), g0[i]], p_ra[i], p_dec[i])
    keep = np.array(keep)
    plt.figure()
    plt.scatter(keep[:,0], keep[:,1], s=1)
    lims = plt.ylim()
    plt.ylim(lims[1], lims[0])
    plt.ylabel(r"$g_0$")
    plt.xlabel(r"$(g-r)_0$")
    plt.show()
    plt.close()

def split_by_plate():
    # dered_g,dered_r,l,b,ELODIERVFINAL,ELODIERVFINALERR,plate,p_ra,p_dec
    data = np.loadtxt("/home/newbym2/Dropbox/Research/sgrLetter/sgr_spec_all.csv", delimiter=",")
    mean, stdev, nsig = -118.233333, 30.222222, 1.0
    plates, pos = [], []
    for i in range(len(data[:,0])):
        if plates.count(data[i,6]) > 0:  continue
        plates.append(data[i,6])
        pos.append([data[i,7], data[i,8]])
    out = []
    for i, plate in enumerate(plates):
        # get stars for plate
        stars = []
        for j in range(len(data[:,0])):
            if data[j,6] == plate:  stars.append(data[j,:])
        if stars == []:  continue
        # apply cuts
        keep = []
        for star in stars:
            if star[0] < 20.0:  continue
            if star[0] > 20.5:  continue  #using stricter cut
            X,Y,Z, lsgr, bsgr, r_sgr = ac.lb2sgr(star[2], star[3], ac.getr(star[0]))
            if lsgr < 190.0:  continue
            if lsgr > 250.0:  continue
            if abs(bsgr) > 10.0:  continue
            vgsr = ac.rv_to_vgsr(star[4],star[2],star[3])
            keep.append([lsgr, bsgr, vgsr])
        if keep == []:  continue
        # vgsr selection
        yes, no = 0, 0
        for k in keep:
            if k[2] < (mean - stdev*nsig):  no += 1;  continue
            if k[2] > (mean + stdev*nsig):  no += 1;  continue
            yes += 1
        print "Plate {0}:  In {1} of Sgr versus out:  {2}, {3}".format(
            plate, nsig, yes, no)
        # compile plate info positions
        out.append([plate, pos[i][0], pos[i][1], yes, no])
    savename = "/home/newbym2/Dropbox/Research/sgrLetter/plate_stars.csv"
    np.savetxt(savename, sc.array(out), delimiter=",")
        

def photo_spec_analysis():
    # l, b, r;  plate #, p_ra, p_dec, in Sgr, out Sgr
    path="/home/newbym2/Desktop/starfiles"
    wd="/home/newbym2/Dropbox/Research/sgrLetter/"
    files = glob.glob(path+"/stars*")
    plates = np.loadtxt(wd+"plate_stars.csv", delimiter=",")
    radius = 1.49
    for i, plate in enumerate(plates[:,0]):
        print "starting plate", str(plate).split(".")[0], plates[i,3], plates[i,4]
        p_ra, p_dec = plates[i,1], plates[i,2]
        pl, pb = ac.EqTolb(p_ra, p_dec)
        print pl, pb
        stars = []
        for f in files:
            data = open(f, "r") #np.loadtxt(f) #, delimiter="\t")
            print "    starting file", f
            wedge = int(f.split("-")[1].split(".")[0])
            skip = 0
            for line in data:
                if skip < 1:  skip =+ 1;  continue
                if line.strip()=="":  continue
                temp = line.split()
                ll, bb, g0 = float(temp[0]), float(temp[1]), ac.getg(float(temp[2]))                
                if g0 < 20.0:  continue
                if g0 > 20.5:  continue
                del_l, del_b = (ll-pl), (bb-pb)
                dd = ma.sqrt( (del_l*del_l) + (del_b*del_b) )
                if dd > radius: continue
                if ac.SDSS_primary(ll,bb,wedge,fmt="lb",low=9,high=27)==0:  continue
                stars.append([ll,bb,g0])
            data.close()
        svname = wd+"plate_"+str(plate).split(".")[0]+"_stars.csv"
        np.savetxt(svname, sc.array(stars), delimiter=",")
        print "saved {0}".format(svname)

def compile_tables():
    # cycle through files using a list of destinations and names
    rnames = ["smarter_slide_back"]
    #rnames = ["smart_sliding_back"]
    """rnames = ["virgo_B-slice_all",
              "back_m2b160_all",
              "flatsub_100_all",
              "cutoff_R30_all",
              "virgo_B-slice_spec",
              "back_m2b160_spec",
              "flatsub_100_spec",
              "back_sliding"
              ]"""
    wd="/home/newbym2/Dropbox/Research/sgrLetter/"
    folder = ["fits_smarter_slide/"]
    #folder = ["fits_smart_slide/"]
    """folder = ["no_virgo_B_slice_smallbins/",
              "backfunc_smallbins/",
              "hist_fits_flatsub100/",
              "hists_r30_cutoff/",
              "no_virgo_B_slice_spec_cut/",
              "backfunc_spec_cut/",
              "hist_fits_spec_cut_flatsub30/",
              "fits_sliding_line"
              ]"""
    subfolder = ["double_gaussian/", "quad_gaussian/"]
    fname = "AA_hist_fits.txt"
    outf = "fit_results/"
    for i, f in enumerate(folder):
        for sf in subfolder:
            run = rnames[i]+"_"+sf.split("_")[0]
            results = open(wd+f+sf+fname, "r")
            table, entry = [], 0
            for line in results:
                if line.strip()=="":  continue
                if (len(line.strip()) < 10):
                    if entry != 0:  table.append(entry) 
                    entry = line.strip().strip(":")
                if line[4:8] == "MCMC":
                    temp = line.split("[")[1].split("]")[0]
                    entry = entry+", "+temp
                if line[4:9] == "Value":
                    temp = line.split("[")[1].split("]")[0]
                    entry = entry+", "+temp
            table.append(entry)
            results.close
            outfile = open(wd+outf+run+".txt", "w")
            for thing in table:
                outfile.write(thing+"\n")
            outfile.close
            #print table #.sort()

def plot_ganged_hists():
    # Make 4x3 blocks of histgram fits;  only have to change next line
    fname = "smarter_slide_back_quad.txt"
    wd = "/home/newbym2/Dropbox/Research/sgrLetter/fit_results/"
    r = np.loadtxt(wd+fname, delimiter=",")
    results = r[np.lexsort( (r[:,0], r[:,0]) )]
    path = "/home/newbym2/Dropbox/Research/sgrLetter/hist_fits_smallbins/"
    #files = glob.glob(path+"/*.out")
    plate_info = np.loadtxt("/home/newbym2/Dropbox/Research/sgrLetter/plate_data.csv", delimiter=",")
    # range of plots
    b_range = (-30, 30)
    N_range = (0, 450)
    # subplot organization
    px, py = 3, 4
    #p1 = [9,5,1,10,11,12,2,3,4,6,7,8]
    # initialize plot info
    if fname[-8:-4] == "quad":  function = func.quad_fat_gauss_line
    else:  function = func.double_gauss_line
    # Panes which skip labelling axes
    nox = [1,2,3,4,5,6,7,8,9]
    noy = [2,3,5,6,8,9,11,12]
    x = np.arange(b_range[0], b_range[1], 0.1)
    # get this bastard going
    for i in range(len(results[:,0])):
        if i == 0:
            fig1 = plt.figure(num=1, figsize=(12,9), dpi=100)
            plt.subplots_adjust(hspace=0.001, wspace=0.001)
            oo = 0  #offset
        if i == 12:  
            fig2 = plt.figure(num=2, figsize=(12,9), dpi=100)
            plt.subplots_adjust(hspace=0.001, wspace=0.001)
            oo = -12
        if i == 24:
            fig3 = plt.figure(num=3, figsize=(12,9), dpi=100)
            plt.subplots_adjust(hspace=0.001, wspace=0.001)
            oo = -24
        if i == 36:
            fig4 = plt.figure(num=4, figsize=(8,4.5), dpi=100)
            plt.subplots_adjust(hspace=0.001, wspace=0.001)
            px, py = 2, 2
            nox, noy = [1,2], [2,4]
            oo = -36
        Lslice = results[i,0]
        # spec stuff
        holder = []
        for j in range(len(plate_info[:,0])):
            if abs(plate_info[j,0]-Lslice) > 1.25:  continue
            holder.append([plate_info[j,0], plate_info[j,1], plate_info[j,2], plate_info[j,3]])
        # get histogram
        hist=np.loadtxt(path+str(results[i,0])+".out")
        #Make plot
        sp = plt.subplot(py,px,(i+1+oo))
        sp.bar(hist[:,0], hist[:,1], 0.5, align='center', color='w', zorder=1)
        if fname[-8:-4] == "quad":
            plt.plot(x, function(x, results[i,1:13]), 'k-', zorder=2)
        else: plt.plot(x, function(x, results[i,1:9]), 'k-', zorder=2)
        plt.plot(x, func.gaussian_function(x, results[i,1:4]), 'k:', zorder=2)
        plt.plot(x, func.gaussian_function(x, results[i,4:7]), 'k:', zorder=2)
        if fname[-8:-4] == "quad":
            plt.plot(x, func.gaussian_function(x, [results[i,7], results[i,2], results[i,8]]), 'k:', zorder=2)
            plt.plot(x, func.gaussian_function(x, [results[i,9], results[i,5], results[i,10]]), 'k:', zorder=2)
        # spec stuff
        if holder != []:
            hold=np.array(holder)
            if fname[-8:-4] == "quad":  norm = function(hold[:,1], results[i,1:13])
            else:  norm = function(hold[:,1], results[i,1:9])
            plt.scatter(hold[:,1], norm*hold[:,2], c="red", marker="s", zorder=5)
            #plt.errorbar(hold[:,1], norm*hold[:,2], yerr= norm*(np.sqrt(1+(hold[:,2]*hold[:,3]))+1), 
            #    marker=None, ls=" ", zorder=4, ecolor="red")
            plt.scatter(hold[:,1], norm, c="green", marker="s", zorder=4)
            #plt.errorbar(hold[:,1], norm, yerr= norm*(np.sqrt(1+((1-hold[:,2])*hold[:,3]))+1), 
            #    marker=None, ls=" ", zorder=3, ecolor="green")
            for k in range(len(hold[:,0])):
                plt.text(hold[k,1], 400, str(int(hold[k,3])), fontsize=8)
        #plt.title(name, fontsize=8 )
        plt.text(-25, 375, r"$\Lambda$="+str(Lslice), fontsize=12)
        plt.xlim(b_range[0], b_range[1])
        plt.ylim(N_range[0], N_range[1])
        # i+1 corresponds to plot number starting from 1
        if (i+1+oo) in nox:  plt.setp(sp.get_xticklabels(), visible=False)
        else:  
            plt.xlabel("B")
            plt.xticks(np.arange(-30.0, 30.0, 5.0), ["-30","","-20","","-10","","0","","10","","20","",""])
        if (i+1+oo) in noy:  plt.setp(sp.get_yticklabels(), visible=False)
        else:  
            plt.ylabel("N")
            plt.yticks(np.arange(0.0, 450.0, 50.0), ["","","100","","200","","300", "","400"])
    plt.show()
    plt.close()

def process_spec():
    pp = np.loadtxt("/home/newbym2/Dropbox/Research/sgrLetter/plate_stars.csv", delimiter=",")
    pl,pb = ac.EqTolb(pp[:,1], pp[:,2])
    pd = 20.0
    pX,pY,pZ, plsgr, pbsgr, pr_sgr = ac.lb2sgr(pl, pb, pd)
    norm = pp[:,3] + pp[:,4]
    out = []
    for i in range(len(pp[:,0])):
        out.append([plsgr[i], pbsgr[i], pp[i,4]/norm[i], norm[i]])
    np.savetxt("/home/newbym2/Dropbox/Research/sgrLetter/plate_data.csv", np.array(out), delimiter=",")

def spec_area():
    """ Returns the ratio of in-Sgr stars to non-Sgr stars 
        within 1-sigma of the Sgr mean, given fits to data within +/-10.0 Beta """
    fit = [3.1707598689912997, -113.54662564418646, 29.872019523468833, 3.539017234542301, 0.0, 120.0]
    x0, xf = (-118.2-30.2), (-118.2+30.2)
    sgr = func.integrate_gaussian(x0, xf, mu=-118.2, sig=30.2, amp=(fit[0]*fit[0]) )
    sgrlost = func.integrate_gaussian(x0-60.4, xf+60.4, mu=-118.2, sig=30.2, amp=(fit[0]*fit[0]) ) - sgr
    notsgr = func.integrate_gaussian(x0, xf, mu=0.0, sig=120.0, amp=(fit[3]*fit[3]) )
    #notsgr = func.integrate_gaussian((-500.0), (500.0), mu=0.0, sig=120.0, amp=(fit[3]*fit[3]) )
    print sgr, notsgr, sgr/notsgr
    print sgrlost, sgrlost/notsgr
    """ 367.393783337 329.186581635 1.1160654894
        169.309797526 0.514327761127"""

def make_heatmap():
    """ produces a heatmap showing probability of star detection for each pixel"""
    fname = "smart_sliding_back_double.txt"
    mask = "Rhist_sgr.csv"
    wd = "/home/newbym2/Dropbox/Research/sgrLetter/fit_results/"
    r = np.loadtxt(wd+fname, delimiter=",")
    fitdata = r[np.lexsort( (r[:,0], r[:,0]) )]
    mask = np.loadtxt(mask_data, delimiter=",")
    lmin, lmax, lstep = 200.0, 300.0, 0.5
    bmin, bmax, bstep = -40.0, 30.0, 0.5
    lbinsize, bbinsize = 2.5, 0.5

def count_stars_in_rcut():
    rmin, rmax = ac.getr(16.0), ac.getr(23.5)
    print rmin, rmax
    path="/home/newbym2/Desktop/starfiles"
    files = glob.glob(path+"/stars*")
    total = 0
    for f in files:
        wedge = int(f.split("-")[1].split(".")[0])
        count = 0
        data = np.loadtxt(f, skiprows=1)
        for i in range(len(data[:,0])):
            if data[i,2] < rmin:  continue
            if data[i,2] > rmax:  continue
            if ac.SDSS_primary(data[i,0],data[i,1],wedge,fmt="lb",low=9,high=27)==0:  continue
            count = count + 1
            total = total + 1
        print count, len(data[:,2]), np.ma.min(data[:,2]), np.ma.max(data[:,2])
    print "Final Count:", total
    

def count_stars():
    fname = "smart_sliding_back_quad.txt"
    wd = "/home/newbym2/Dropbox/Research/sgrLetter/fit_results/"
    r = np.loadtxt(wd+fname, delimiter=",")
    data = r[np.lexsort( (r[:,0], r[:,0]) )]
    path = "/home/newbym2/Dropbox/Research/sgrLetter/hist_fits_smallbins/"
    sgr1, sgr2, virgo, back = 0, 0, 0, 0
    if fname[-8:-4] == "quad":  sgr1b, sgr2b = 0, 0
    dt, dL = 0.5, 2.5
    all_stars = 0
    for i in range(len(data[:,0])):
        hist = np.loadtxt(path+str(data[i,0])+".out")
        all_stars = all_stars + sum(hist[:,1])*dL
        #print np.sum(hist[:,1])
        #get range of data
        j,x = -1, 0  
        while x < 100:  x = hist[j,1]; j -= 1
        high = hist[j,0]
        j,x = 0,0  
        while x < 100:  x = hist[j,1]; j += 1
        low = hist[j,0]
        #Integrate to get total stars
        if fname[-8:-4] == "quad":  aa, bb = data[i,11], data[i,12]
        else:  aa, bb = data[i,7], data[i,8]
        if aa < 0:
            newback = (aa*low + bb)*(high-low)/dt
            newvirgo = 0
        else:
            newback = (aa*low + bb)*(high-low)/dt
            newvirgo = (aa*high - aa*low)*(high - low)/(2.0*dt)
        newsgr1 = func.integrate_gaussian(low, high, mu=data[i,2], sig=data[i,3], 
            amp=(data[i,1]*data[i,1]) )/dt
        newsgr2 = func.integrate_gaussian(low, high, mu=data[i,5], sig=data[i,6], 
            amp=(data[i,4]*data[i,4]) )/dt
        if fname[-8:-4] == "quad":
            newsgr1b = func.integrate_gaussian(low, high, mu=data[i,2], sig=data[i,8], 
                amp=(data[i,7]*data[i,7]) )/dt
            newsgr2b = func.integrate_gaussian(low, high, mu=data[i,5], sig=data[i,10], 
                amp=(data[i,9]*data[i,9]) )/dt
        back = back + abs(newback)*dL
        virgo = virgo + abs(newvirgo)*dL
        sgr1 = sgr1 + abs(newsgr1)*dL
        sgr2 = sgr2 + abs(newsgr2)*dL
        if fname[-8:-4] == "quad":
            sgr1b = sgr1b + abs(newsgr1b)*dL
            sgr2b = sgr2b + abs(newsgr2b)*dL
    if fname[-8:-4] == "quad":
        total = (sgr1+sgr1b+sgr2+sgr2b+virgo+back)
        print "Sgr1a: {0:8.1f},  {1:.3f}".format(sgr1, sgr1/total)
        print "Sgr1b: {0:8.1f},  {1:.3f}".format(sgr1b, sgr1b/total)
        print "Sgr2a: {0:8.1f},  {1:.3f}".format(sgr2, sgr2/total)
        print "Sgr2b: {0:8.1f},  {1:.3f}".format(sgr2b, sgr2b/total)
        print "Virgo: {0:8.1f},  {1:.3f}".format(virgo, virgo/total)
        print "Back : {0:8.1f},  {1:.3f}".format(back, back/total)
        print "Total: {0:8.1f}".format(total)
    else:
        total = (sgr1+sgr2+virgo+back)
        print "Sgr1 : {0:8.1f},  {1:.3f}".format(sgr1, sgr1/total)
        print "Sgr2 : {0:8.1f},  {1:.3f}".format(sgr2, sgr2/total)
        print "Virgo: {0:8.1f},  {1:.3f}".format(virgo, virgo/total)
        print "Back : {0:8.1f},  {1:.3f}".format(back, back/total)
        print "Total: {0:8.1f}".format(total)
    print "All Stars: {0}".format(int(all_stars))
    

def make_data_tables():
    fname = "smarter_slide_back_quad.txt"
    wd = "/home/newbym2/Dropbox/Research/sgrLetter/fit_results/"
    r = np.loadtxt(wd+fname, delimiter=",")
    data = r[np.lexsort( (r[:,0], r[:,0]) )]
    #set offset index for errors
    if fname[-8:-4] == "quad":  
        de = 12
        print "slice\t A1\t mu1\t sig1\t A2\t mu2\t sig2\t A1b\t sig1b\t A2b\t sig2b\t aa\t bb"
    else:  
        de = 8
        print "slice\t A1\t mu1\t sig1\t A2\t mu2\t sig2\t aa\t bb"
    #print data[0,:]
    pm = u"\u00B1"  #unicode plus-minus symbol
    for i in range(len(data[:,0])):
        out = ""
        out = out+str(data[i,0])+" "
        for j in range(1, de+1):
            out = out+"{0:6.1f}".format(data[i,j])+pm+"{0:.2f}".format(data[i,(j+de)])
        print out
        
def sgr_plot3D():
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    fname = "smarter_slide_back_quad.txt"
    wd = "/home/newbym2/Dropbox/Research/sgrLetter/fit_results/"
    r = np.loadtxt(wd+fname, delimiter=",")
    data = r[np.lexsort( (r[:,0], r[:,0]) )]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    z = np.arange(200.0, 300.0, 2.5)    
    virgox, virgoy = [], []
    for i in range(len(data[:,0])):
        xs = np.arange(-30.0, 30.0, 0.1)
        if fname[-8:-4] == "quad":  
            ys = func.quad_fat_gauss_line(xs, data[i,1:13])
        else:  
            ys = func.double_gauss_line(xs, data[i,1:9])
        zz = z[i] #(z[i]-200.0)*10.0
        ax.plot(xs, ys, zs=zz, zdir='y', alpha=0.8, color="black")
        virgox.append(zz)
        virgoy.append(ys[-1])
    """# build contour plot
    path = "/home/newbym2/Dropbox/Research/sgrLetter/hist_fits_smallbins"
    files = glob.glob(path+"/*.out")
    X, Y, Z = [], [], []
    for f in files:
        name = f[-9:-4]
        zslice = float(name)
        hist = np.loadtxt(f)
        tempy, tempz = [], []
        X.append(hist[:,0])
        for j in range(len(hist[:,0]) ):
            tempy.append(zslice)
            tempz.append(hist[j,1])
        Y.append(tempy)
        Z.append(tempz)
    X, Y, Z = np.array(X), np.array(Y), np.array(Z)
    Z = np.sqrt(Z)
    con1 = ax.contourf(X, Y, Z, zdir='z', offset=0.0, cmap=pp.spectral_wb) """
    #    ax.bar(hist[:,0], hist[:,1], zs=float(name), zdir='y', alpha=0.8, color="white")
    #vout = np.transpose(np.vstack((np.array(virgox), np.array(virgoy) ) ) )
    #print vout, vout.shape
    #np.savetxt("virgo_edge_points.txt", vout, fmt='%.6f') 
    ax.plot(virgox, virgoy, zs=30.0, zdir='x', alpha=0.8, color="red")
    ax.plot(virgox, func.gauss_plus_floor(virgox, [13.0, 267., 20.5, 137.]), zs=30.0, zdir='x', alpha=0.8, color="cyan")
    ax.set_xlabel('B')
    ax.set_ylabel(r'$\Lambda$')
    ax.set_zlabel('N')
    ax.set_xlim(-30.0, 30.0)
    ax.set_ylim(200.0, 300.0)
    ax.set_zlim(0.0, 600.0)
    plt.show()

def tomography():
    cutstep = 0.5
    gcuts = np.arange(16.0, 24.0, cutstep)
    for i in range(len(gcuts)):
        rmin, rmax = ac.getr(gcuts[i]), ac.getr(gcuts[i]+cutstep)
        outname = "SDSSN_gslice_{0:.1f}_{1:.1f}.png".format(rmin, rmax)
        make_total_plot(RGB=0, rcut=(rmin, rmax), outfile=outname )
    print "DONE"

def crotus_cut(path="/home/newbym2/Desktop/starfiles", cdata=None):
    if cdata == None:
        files = glob.glob(path+"/stars*")
        data=[]
        r_cut_low, r_cut_high = ac.getr(16.0), ac.getr(23.5) #30.0, 45.0
        pb = pr.Progressbar(steps=len(files), prefix="Loading Stars:", suffix=None,
            symbol="#", active="=", brackets="[]", percent=True, size=40)
        for f in files:
            wedge = int(f.split("-")[1].split(".")[0])
            sdata = np.loadtxt(f, skiprows=1)
            for i in range(sdata.shape[0]):
                if (sdata[i,2] < r_cut_low) or (sdata[i,2] > r_cut_high):  continue
                if ac.SDSS_primary(sdata[i,0],sdata[i,1],wedge,fmt="lb",low=9,high=27)==0:  continue
                lam, bet = (ac.lb2sgr(sdata[i,0], sdata[i,1], 10.0))[3:5]
                if (lam <230.0) or (lam > 255.0):  continue
                data.append([sdata[i,0], sdata[i,1], sdata[i,2], lam, bet])
            pb.updatebar(float(files.index(f)+1)/float(len(files)) )
        pb.endbar()
        data = np.array(data)
        np.savetxt("crotus_data.txt", data, fmt='%.6f')
    else:  
        data = []
        rdata = np.loadtxt(cdata)
        rlim0 = 28.5  #ac.getr(22.5)
        rlim1 = 45.0
        for i in range(rdata.shape[0]):
            if rdata[i,2] < rlim0:  continue
            if rdata[i,2] > rlim1:  continue
            data.append(rdata[i,:])
        data = np.array(data)
    #Now do analysis
    sky = pp.HistMaker(data[:,3], data[:,4], xsize=0.25, ysize=0.25, 
        xarea=(230.0, 255.0), yarea=(-10.0, 15.0))
    sky.varea = (0.0, 60.0)
    sky.cmap="color"
    #sky.scale = 'sqrt'
    sky.yflip = 1
    pp.PlotHist(sky, "crotus_cut.png")
    #sky.savehist("streamgen_bifExtra.csv")
    print "### - Done"
          
def lam_wedges():
    sys.path.insert(0, '../milkyway-tools')
    import sdss_visualizers as sdss
    path = "/home/newbym2/Desktop/lam_starfiles/"
    files = glob.glob(path+"/stars*")
    for f in files:
        # get lambdas from naming scheme
        wedge = float(f[-6:-4])
        name = "lam-"+str(200.0 + (wedge*2.5) )[:6]
        stars = np.loadtxt(f)
        print "Loaded ", f
        for i in range(stars.shape[0]):
            stars[i,0], stars[i,1] = (ac.lb2sgr(stars[i,0], stars[i,1], 30.0))[3:5]
        sdss.plot_stripe_mur(stars, wedge, outname=name, mag=1, scale=1, color=1,
                    mu_lim=None, r_lim=(0.0, 50.0), vm=10.0, nu_flatten=0, bar=1, raw_coords=1)
        print "Finished file {0} of {1}:".format(files.index(f)+1, len(files) )
    print "### - Done"
            
if __name__ == "__main__":
    #shift_sgr(filein="streamgen_sgrfidprim.txt", fileout="stream_shiftfid.txt")
    #make_sim_stream_plot(filein="streamgen_sfp_bigish.txt", RGB=1) #, imfile="sgr_new.png")
    #make_total_plot(RGB=0)
    #make_total_plot(RGB=0, rcut=(ac.getr(20.0), ac.getr(20.5) ) )
    #make_diff_hist()
    #get_bif()
    #get_sgr_curves()
    #batch_shift()
    #RGB_from_files(mask_data="Rhist_sgr.csv", imfile="new.png")
    #RGB_from_files(mask_data="Rhist_sgr.csv", imfile=None, fitfile=True)
    plot_profiles()    
    #plot_profiles(savewedge=True)
    #proj_test()
    #sgr_rv_skyplot()
    #sgr_rv_cut()
    #split_by_plate()
    #photo_spec_analysis()
    #compile_tables()
    #process_spec()
    #plot_ganged_hists()
    #spec_area()
    #make_data_tables()
    #count_stars()
    #count_stars_in_rcut()
    #sgr_plot3D()
    #tomography()
    #crotus_cut(cdata="crotus_data.txt")
    #lam_wedges()
    
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
    
"""
