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
import progress as pr
import glob

def make_sim_stream_plot():
    """ Makes the plot for the simulated streams 
        /home/newbym2/Dropbox/Research/sgrnorth_paper/sgr_separated_stars_MRT.txt"""
    filename="streamgen_bifm10_sim.txt"
    file2="streamgen_sgr_sim.txt"
    data = np.loadtxt(filename)
    for i in range(len(data[:,0])):
        data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
    data2 = np.loadtxt(file2)
    for i in range(len(data2[:,0])):
        data2[i,0], data2[i,1] = ac.lbToEq(data2[i,0], data2[i,1])
    sky = pp.HistMaker(np.concatenate([data[:,0],data2[:,0]]), np.concatenate([data[:,1],data2[:,1]]), 
        xsize=0.5, ysize=0.5, xarea=(120.0, 250.0), yarea=(-10.0, 50.0))
    sky.varea = (0.0,200.0)
    sky.H = sky.H + np.random.normal(45.0, 15.0, sky.H.shape)
    #for i in range(sky.H.shape[0]):
    #    if i < 14:   sky.H[i,:] = sky.H[i,:]*0.0; continue
    #    sky.H[i,:] = sky.H[i,:] + 45.0 - (0.5/1.0)*i
    pp.PlotHist(sky, "simtotbifm10_radec.png")
    sky.savehist("simtotbifm10_radec.csv")
    print "Ended Successfully"

def make_total_plot(path="/home/newbym2/Desktop/starfiles"):
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
            if test_primary(temp[0],temp[1],wedge,low=9,high=27)==0:  continue
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
        data[i,0], data[i,1] = ac.lbToEq(data[i,0], data[i,1])
        if count % 100 == 0:  pb2.updatebar(float(count)/nStars)
    pb2.endbar()
    allsky = pp.HistMaker(data[:,0], data[:,1], xsize=0.5, ysize=0.5, 
        xarea=(120.0, 250.0), yarea=(-10.0, 50.0))
    #allsky.scale = 'sqrt'
    allsky.varea = (0.0,200.0)
    pp.PlotHist(allsky, "sgrall_GC.png")
    allsky.savehist("SDSSnorthGC.csv")
    print "Ended Successfully"

def make_diff_hist():
    """ Subtracts the simulated streams from the real data """
    data = np.loadtxt("SDSSnorth.csv", delimiter=",")
    #sim = np.loadtxt("simtotbif50m30_radec.csv", delimiter=",")
    sim = np.loadtxt("simtotbifm2_radec.csv", delimiter=",")
    #new = data-sim
    for i in range(sim.shape[0]):
        for j in range(sim.shape[1]):
            if data[i,j] < 1.0:  sim[i,j] = 0.0
    sky = pp.HistMaker([1,2], [1,2], xsize=0.5, ysize=0.5, 
        xarea=(120.0, 250.0), yarea=(-10.0, 50.0))
    sky.H = sim
    sky.varea = (0.0,200.0)
    #sky.plot()
    pp.PlotHist(sky, "simbifm_zoomed.png")
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
        
    
def test_primary(l,b,wedge,low=9,high=23):
    """ Tests to see if a star is primary for its wedge number, by testing its
        nu against adjecent stripes """
    mu0, nu0 = ac.lb2GC(l,b,wedge)
    if wedge > low:
        mu1, nu1 = ac.lb2GC(l,b,wedge-1)
        if abs(nu1) < abs(nu0):  return 0
    if wedge < high:
        mu2, nu2 = ac.lb2GC(l,b,wedge+1)
        if abs(nu2) < abs(nu0):  return 0
    return 1

if __name__ == "__main__":
    #make_sim_stream_plot()
    #make_total_plot()
    #make_diff_hist()
    get_bif()

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
