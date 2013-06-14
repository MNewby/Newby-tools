import sys
sys.path.insert(0, '../utilities')
import math as ma
import numpy as np
import scipy as sc
import files as fi
import fileinput
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
#from mpl_toolkits.mplot3d import Axes3D
import time
import astro_coordinates as coor
import os.path
import sys
import glob
import test_wedge_generator as twg

'''Input: l,b,g_0
l limits: 20.4835409501, 205.337883034
b limits: -85.9538084201, -30.0000008077
g limits: 16, 23
Matthew Newby, August 30, 2011'''

Sgr_lbr = [5.6, -14.2, 24.0]
Sgr = coor.lbr2xyz(Sgr_lbr[0], Sgr_lbr[1], Sgr_lbr[2])

class stripe:
    """ Do not include extension in file name!"""
    def __init__(self, filename=""):
        self.filename = filename
        self.data = []
        self.size = 0
        self.file_size = 0
        while os.path.isfile(filename) == True:
            filename = filename + "_new"
        filename = filename + ".txt"
        file = open(filename, "w")
        file.write("#  l, b, r")
        file.close()
    def add(self, input):
        self.data.append(input)
        self.size = self.size + 1
        if self.size > 9999:
            self.dump()
    def dump(self):
        if len(self.data) > 0:
            fi.append_data(sc.array(self.data), self.filename)
        self.data = []
        self.file_size = self.file_size + self.size
        self.size = 0

deg = 180.0 / ma.pi 
rad = ma.pi / 180.0 

# imshow counts from top!!!

cdict1 = {'red':  ((0.0, 1.0, 1.0),
                   (0.25, 0.75, 0.75),
                   (0.5, 0.5, 0.5),
                   (0.75, 0.25, 0.25),
                   (1.0, 0.0, 0.0)),
         'green': ((0.0, 1.0, 1.0),
                   (0.25, 0.75, 0.75),
                   (0.5, 0.5, 0.5),
                   (0.75, 0.25, 0.25),
                   (1.0, 0.0, 0.0)),
         'blue':  ((0.0, 1.0, 1.0),
                   (0.25, 0.75, 0.75),
                   (0.5, 0.5, 0.5),
                   (0.75, 0.25, 0.25),
                   (1.0, 0.0, 0.0))        }
white_black = LinearSegmentedColormap('grey', cdict1)
spectral_colors = fi.read_data('../utilities/spectral_cm.txt')
spectral_wb = ListedColormap(spectral_colors[:,:3], name="spectral_wb", N=256)

def get_stripe(wedge, filein):
    stripe_data = []
    skymap = sc.zeros((180, 360), float)
    for line in fileinput.input(filein):
        # Get Data
        if line.strip()[0] == '':  continue
        if line.strip()[0] == '#': continue
        holder = line.split(',')
        l,b,g = eval(holder[0].strip()), eval(holder[1].strip()), eval(holder[2].strip())
        ra, dec = coor.lbToEq(l,b)
        mu, nu = coor.EqToGC(ra, dec, wedge)
        if abs(nu) < 1.25:
            stripe_data.append([l,b,g])
            # now build the histogram
            l_i = int(ra)
            b_i = int(dec) + 90
            skymap[-b_i][l_i] = skymap[-b_i][l_i] + 1.0
        if fileinput.filelineno() % 100000 == 0:
            print "Line {0} complete".format(fileinput.filelineno())
    return sc.array(stripe_data), skymap
    
def plot_sky(skymap, wedge, out_form="stripe_"):
    plotout = out_form+str(wedge)+".ps"
    extent = [0.0, 360.0, -90.0, 90.0]
    mu = sc.arange(-90.0, 90.0, 10.0)
    nu = sc.zeros(len(mu))
    ra, dec = coor.GCToEq(mu, nu, wedge)
    plt.figure(1)
    plt.imshow(skymap, extent=extent)
    plt.scatter(ra,dec,s=2,c="r")
    plt.savefig(plotout, papertype='letter')
    #plt.show()
    
def all_stripes(stripes, file_in, out_form="stripe_"):
    t1 = time.time()
    for stripe in stripes:
        fileout = out_form+str(stripe)+".txt"
        out, map = get_stripe(stripe, file_in)
        fi.write_data(out, fileout)
        plot_sky(map, stripe)
        t2 = time.time()
        print "Time elapsed, stripe {0}: {1} minutes".format(stripe, ((t2-t1)/60.0))

def plot_stripe_lb(data, bin_size=0.5):
    l_min, l_max = np.ma.min(data[:,0]), np.ma.max(data[:,0])
    b_min, b_max = np.ma.min(data[:,1]), np.ma.max(data[:,1])
    l_bins, b_bins = (int((l_max-l_min)/bin_size)+1), (int((b_max-b_min)/bin_size)+1)
    H, x, y = np.histogram2d(data[:,0], data[:,1], [l_bins, b_bins])
    extent = [b_min, b_max, l_min, l_max]
    plt.imshow(H, extent=extent, interpolation='nearest')
    plt.colorbar()
    plt.show()
    
def plot_stripe_3D(data):
    l, b, r = data[:,0], data[:,1], coor.getr(data[:,2])
    x,y,z = coor.lbr2xyz(l,b,r)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter([-8.5], [0.0], [0.0], c='red')
    ax.scatter(x,y,z, c='yellow')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    
""" --------------- SGR Plot ------------------------------"""
def sgr_plot(folder="./sep_lbr"):
    files = glob.glob(folder+"/*.txt")
    plt.figure(1)
    for file in files:
        data = fi.read_data(file)
        Xs, Ys, Zs, lam, beta, r = coor.lb2sgr(data[:,0], data[:,1], data[:,2])
        plt.scatter(lam, beta, s=1, c="black")
    plt.show()
    plt.close('all')

""" --------------- Side-View Sgr plot (XZ) --------------- """
def sgr_xz_plot(folder="./sep_lbr/", bin_size=1.0, outfile=None, infile=None,
                primary=1, color=1, scale=1):
    x_min, x_max, x_size = -50.0, 50.0, bin_size
    z_min, z_max, z_size = 0.0, 100.0, bin_size
    x_bins, z_bins = int((x_max-x_min)/x_size)+1, int((z_max-z_min)/z_size)+1
    if infile==None:
        H = sc.zeros((z_bins,x_bins), float)
        files = glob.glob(folder+"/*.txt")
        for file in files:
            data = fi.read_data(file)
            wedge = int(file[-6:-4])
            for i in range(len(data[:,0])):
                if primary==1:
                    if test_primary(data[i,0],data[i,1],wedge,low=9,high=23)==0:  continue
                x,y,z = coor.lbr2xyz(data[i,0], data[i,1], data[i,2])
                if x < x_min or x > x_max:  continue
                if z < z_min or z > z_max:  continue        #FLAG
                xi, zi = int((x-x_min)/x_size), int((z-z_min)/z_size) #FLAG
                H[zi,xi] = H[zi,xi] + 1.0
        if outfile!=None:  fi.write_data(H, outfile)
    else:  H = fi.read_data(infile)
    if scale==1:  H = sc.sqrt(H)
    #Make the plot
    plt.figure(1)
    sp = plt.subplot(111)
    if color==1:  cmap = spectral_wb
    else:  cmap = binary
    plt.imshow(H, interpolation='nearest', cmap=cmap) #vmin=0.0, vmax=15.0)
    bar = plt.colorbar(orientation='vertical')
    if scale == 1:
        bar_ticks = sc.arange(0.0, 17.0, 4.0)
        bar_labels = []
        for i in range(len(bar_ticks)):   
            bar_labels.append(str(int(bar_ticks[i]**2)))
        bar.set_ticks(bar_ticks)
        bar.set_ticklabels(bar_labels, update_ticks=True)
    Earth = [-8.5, 0.0, 0.0]
    plt.scatter( (Earth[0]-x_min)/x_size, (Earth[2]-z_min)/z_size, c='yellow',s=5 )
    # Clean up
    y_lim = plt.ylim()
    plt.ylim(y_lim[1], y_lim[0])
    #plt.xlim(-25.0, 270.0)
    #plt.ylim(-25.0, 270.0)
    #plt.setp(sp.get_xticklabels(), visible=False)
    #plt.setp(sp.get_yticklabels(), visible=False)
    #plt.setp(sp.get_xticklines(),  visible=False)
    #plt.setp(sp.get_yticklines(),  visible=False)
    # Show plot
    plt.show()
    plt.close('all')


def plot_stripe_mug(data, wedge, name="out", mag=0):
    """ SPREAD OUT MORE!!! """
    #l, b = data[:,0], data[:,1]
    if mag==1:  g = data[:,2]
    else:       g = coor.getg(data[:,2])
    ra, dec = coor.lbToEq(data[:,0],data[:,1])
    mu, nu = coor.EqToGC(ra, dec, wedge)
    x, y  = g*sc.cos(mu*rad), g*sc.sin(mu*rad)
    x_bins = int((np.ma.max(x)-np.ma.min(x))/0.5) + 1
    y_bins = int((np.ma.max(y)-np.ma.min(y))/0.5) + 1
    H, x, y = np.histogram2d(x, y, [x_bins, y_bins])
    extent = [np.ma.min(y), np.ma.max(y), np.ma.max(x), np.ma.min(x)]
    plt.figure(1)
    sp = plt.subplot(111)
    #plt.polar([0.0, np.ma.min(mu)*rad, np.ma.max(mu)*rad, 0.0], [0.0, 60.0, 60.0, 0.0])
    plt.imshow(H, extent=extent, interpolation='nearest')
    plt.colorbar(orientation='horizontal')
    #plt.plot(x1,y1,'r-')
    plt.setp(sp.get_xticklabels(), visible=False)
    plt.setp(sp.get_yticklabels(), visible=False)
    plt.savefig((name+".ps"), papertype='letter')
    plt.close('all')
    #plt.show()

""" -------------  l,b POLAR PLOT ------------------------------------------ """

def lbpolar_plot(folder, hemi='N', bin_size=1.0, outfile=None, infile=None,
                 color=1, scale=1, primary=1):
    x_size, y_size, b_min = bin_size, bin_size, 30.0
    b_min = 90.0-b_min  # Flips it around so that 90.0 is at the origin
    # Initialize custom histogram
    x_bins, y_bins = int(2.0*b_min/x_size)+1, int(2.0*b_min/y_size)+1
    if infile==None:
        H = sc.zeros((y_bins,x_bins), float)
        # Iterate over files
        files = glob.glob(folder+"/*.txt")
        for file in files:
            data = fi.read_data(file)
            wedge = int(file.split("/")[-1][3:5])  
            for i in range(len(data[:,0])):
                if primary==1:
                    if test_primary(data[i,0],data[i,1],wedge,low=9,high=23)==0:  continue
                xi = int(((90.0-data[i,1])*sc.cos(data[i,0]*rad) + b_min)/x_size)
                yi = int(((90.0-data[i,1])*sc.sin(data[i,0]*rad) + b_min)/y_size)
                if xi < 0 or xi > x_bins-1:  continue
                if yi < 0 or yi > y_bins-1:  continue
                H[yi,xi] = H[yi,xi] + 1.0
        if outfile!=None:  fi.write_data(H, outfile)
    else:  H = fi.read_data(infile, delimiter=None)
    if scale==1:  H = sc.sqrt(H)
    # plot data
    plt.figure(1)
    sp = plt.subplot(111)
    if color==1:  cmap = spectral_wb
    else:  cmap = binary
    plt.imshow(H, interpolation='nearest', cmap=cmap, vmin=0.0, vmax=12.0) #extent=extent, 
    bar = plt.colorbar(orientation='vertical')
    if scale == 1:
        bar_ticks = sc.arange(0.0, 17.0, 4.0)
        bar_labels = []
        for i in range(len(bar_ticks)):   
            bar_labels.append(str(int(bar_ticks[i]**2)))
        bar.set_ticks(bar_ticks)
        bar.set_ticklabels(bar_labels, update_ticks=True)
    # plot axes - l
    for i in [30.0, 45.0, 60.0, 75.0]:
        l = sc.arange(0.0, 361.0, 1.0)
        x = ((90.0-i)/(bin_size))*sc.cos(l*rad)+(b_min/x_size)          # FLAG x_size
        y = ((90.0-i)/(bin_size))*sc.sin(l*rad)+(b_min/y_size)          # FLAG y_size
        if i%2==0:  plt.plot(x,y,'k-')
        else:  plt.plot(x,y,'k:', linewidth=0.5)
        plt.text(x[89]+1.0, y[89]+1.0, r"$b="+str(int(i))+r"^{\circ}$", horizontalalignment='left',
                 verticalalignment='bottom', fontsize=12)
    # plot axes - b
    for i in sc.arange(0.0, 180.0, 15.0):
        b = sc.arange(-1.0*(b_min+5.0), (b_min+5.0)+1.0, 10.0) / bin_size           # FLAG
        x = b*sc.cos(i*rad)+(b_min/x_size)        # FLAG x_size
        y = b*sc.sin(i*rad)+(b_min/y_size)        # FLAG y_size
        if i%90==0:  plt.plot(x,y, 'k-')
        else:  plt.plot(x,y, 'k:')
        if i%30==0:
            if i==0.0:  pre = r"$l="
            else:  pre = r"$"
            x0 = (b[0]-4.0)*sc.cos(i*rad)+(b_min/x_size)          # FLAG x_size
            y0 = (b[0]-4.0)*sc.sin(i*rad)+(b_min/y_size)          # FLAG y_size
            plt.text(x0, y0, r"$"+str(int(i+180.0))+r"^{\circ}$", horizontalalignment='center',
                 verticalalignment='center', fontsize=12, rotation=(i+90.0) )
            x1 = (b[-1]+4.0)*sc.cos(i*rad)+(b_min/x_size)          # FLAG x_size
            y1 = (b[-1]+4.0)*sc.sin(i*rad)+(b_min/y_size)          # FLAG y_size
            plt.text(x1, y1, pre+str(int(i))+r"^{\circ}$", horizontalalignment='center',
                 verticalalignment='center', fontsize=12, rotation=(i-90.0) )
    # Clean up
    #y_lim = plt.ylim()
    #plt.ylim(y_lim[1], y_lim[0])
    plt.xlim(-25.0, 270.0)
    plt.ylim(-25.0, 270.0)
    plt.setp(sp.get_xticklabels(), visible=False)
    plt.setp(sp.get_yticklabels(), visible=False)
    plt.setp(sp.get_xticklines(),  visible=False)
    plt.setp(sp.get_yticklines(),  visible=False)
    # Show plot
    plt.show()
    plt.close('all')
    
def test_primary(l,b,wedge,low=9,high=23):
    """ Tests to see if a star is primary for its wedge number, by testing its
        nu against adjecent stripes """
    mu0, nu0 = coor.lb2GC(l,b,wedge)
    if wedge > low:
        mu1, nu1 = coor.lb2GC(l,b,wedge-1)
        if abs(nu1) < abs(nu0):  return 0
    if wedge < high:
        mu2, nu2 = coor.lb2GC(l,b,wedge+1)
        if abs(nu2) < abs(nu0):  return 0
    return 1

def plot_separation_mur_OLD(datas, wedge, outname=None, mag=0, scale=0, color=1,
                    mu_lim=None, r_lim=None, nu_flatten=0):
    """ Datas is a list of separated stars: [background, stream1, stream2, stream3]
        Use 'None' for a stream if you don't want it plotted. """
    hists, counter = [], 0
    for data in datas:
        if data==None:  continue
        x_size, y_size = 0.5, 0.5
        if mag==1:  r = coor.getr(data[:,2])
        else:       r = data[:,2]
        ra, dec = coor.lbToEq(data[:,0],data[:,1])
        mu, nu = coor.EqToGC(ra, dec, wedge)
        if counter==0:  #Initializes wedge limits for background, if not already done
            if r_lim == None:  r_lim = (np.ma.min(r), np.ma.max(r))
            if mu_lim==None:  mu_lim = (np.ma.max(mu), np.ma.min(mu))
        x, y = r*sc.cos(mu*rad), r*sc.sin(mu*rad)
        x_bins = int((np.ma.max(x)-np.ma.min(x))/x_size) + 1
        y_bins = int((np.ma.max(y)-np.ma.min(y))/y_size) + 1
        H, x, y = np.histogram2d(y, x, [y_bins, x_bins]) #, range=[[-155.0, 110.0],[-190.0, 0.0]])
        if counter==0:
            extent = [np.ma.min(y), np.ma.max(y), np.ma.max(x), np.ma.min(x)]
        if nu_flatten == 1:
            print "!!! DOUBLE-CHECK NU FLATTENING! {0}".format(H.shape)
            for i in range(H.shape[0]):
                for j in range(H.shape[1]):
                    dist = ma.sqrt( ((i-(H.shape[0]/2.0))*y_size)**2 + (j*x_size)**2)
                    H[i,j] = H[i,j] / (dist*2.5)
            H = sc.sqrt(H) #sc.log10(H)
            #vm = 1.5
        elif scale==1:  H=sc.sqrt(H)
        else:  vm=300.0
        hists.append(H)
        counter=1
    """ Setup Figure """
    if nu_flatten==1 or scale==1:  vm=16.0
    else:  vm = 300.0
    if color==1:  cmap = spectral_wb
    else:  cmap = binary
    """ Begin Figure """
    plt.figure(1)
    plt.subplots_adjust(hspace=0.001, wspace=0.001)
    place = [221, 222, 223, 224]
    for i,H in list(enumerate(hists)):
        sp = plt.subplot(place[i], adjustable='box', aspect=1.0)
        plt.imshow(H, extent=extent, interpolation='nearest', vmin=0.0, vmax=vm,
               cmap=cmap)
        if i==3:  #last plot
            bar = plt.colorbar(orientation='vertical')
            if scale == 1:
                bar_ticks = sc.arange(0.0, vm+1.0, 2.0)
                bar_labels = []
                for i in range(len(bar_ticks)):   
                    bar_labels.append(str(int(bar_ticks[i]**2)))
                bar.set_ticks(bar_ticks)
                bar.set_ticklabels(bar_labels, update_ticks=True)
        """ Draw axes - mu """
        mu_axis = gen_mu_axis(mu_lim, r_lim[1])
        plt.plot(mu_axis[0], mu_axis[1], 'k-')
        for i in range(len(mu_axis[2])):
            x0 = [1.0*r_lim[1]*sc.cos(mu_axis[2][i]*rad), 1.1*r_lim[1]*sc.cos(mu_axis[2][i]*rad)]
            y0 = [1.0*r_lim[1]*sc.sin(mu_axis[2][i]*rad), 1.1*r_lim[1]*sc.sin(mu_axis[2][i]*rad)]
            plt.plot(x0,y0, 'k-')
            if mu_axis[3][i] <= 180.0:  rot = -90.0 + mu_axis[3][i]
            else:  rot =  -360.0 + mu_axis[3][i] - 90.0
            x00, y00 = 1.12*r_lim[1]*sc.cos(mu_axis[2][i]*rad), 1.12*r_lim[1]*sc.sin(mu_axis[2][i]*rad)
            plt.text(x00, y00, str(mu_axis[3][i]), rotation=rot, fontsize=10,
                 horizontalalignment='center') #float(mu_axis[3][i]) )
        """ Draw axes - r """
        r_axis = gen_r_axis(mu_lim, r_lim)
        plt.plot(r_axis[0], r_axis[1], 'k-')
        mu1 = sc.arange(mu_lim[0]-6.0, mu_lim[1]+7.0, 1.0)
        for tick in r_axis[2]:
            x1, y1 = tick*sc.cos(mu1*rad), tick*sc.sin(mu1*rad)
            if mu1[-1] < 90.0 or mu1[-1] > 360.0:
                x2, y2 = tick*sc.cos((mu1[-1])*rad), tick*sc.sin((mu1[-1])*rad)
                hl = 'right'
            else:  x2, y2 = tick*sc.cos((mu1[1])*rad), tick*sc.sin((mu1[1])*rad); hl='left'
            plt.plot(x1[:7], y1[:7], "k-")
            plt.plot(x1[-7:], y1[-7:], "k-")
            plt.text(x2, y2, str(int(tick)), rotation=0.0, fontsize=10,
                 horizontalalignment=hl, verticalalignment='bottom')
        # Draw axes - g (TO BE DONE)
        """ Clean up output """
        y_lim = plt.ylim()
        plt.ylim(y_lim[1], y_lim[0])
        plt.setp(sp.get_xticklabels(), visible=False)
        plt.setp(sp.get_yticklabels(), visible=False)
        plt.setp(sp.get_xticklines(),  visible=False)
        plt.setp(sp.get_yticklines(),  visible=False)
    """ Draw Plot """
    if outname == None:  plt.show()
    else:  plt.savefig((outname+".ps"), papertype='letter')
    plt.close('all')

""" ----------------- plot_stripe_mur ----------------------"""

def plot_stripe_mur(data, wedge, outname=None, mag=0, scale=0, color=1,
                    mu_lim=None, r_lim=None, vm=None, nu_flatten=0, bar=1):
    fig = plt.figure(1, frameon=False)
    sp = single_stripe_mur(data, wedge, mag, scale, color, 111, mu_lim, r_lim, vm, nu_flatten, bar)
    """ Draw Plot """
    if outname == None:  plt.show()
    else:  plt.savefig((outname+".ps"), papertype='letter')
    plt.close('all')


def plot_separation_mur(wedge, data0, data1=None, data2=None, data3=None, 
                        outname=None, mag=0, scale=0, color=1, mu_lim=None, 
                        r_lim=None, vm=None, nu_flatten=0, bar=1):
    """ """
    fig = plt.figure(1, frameon=False)
    fig.subplots_adjust(hspace=0.001, wspace=0.001)
    #place = [221, 222, 223, 224]
    sp0 = single_stripe_mur(data0, wedge, mag, scale, color, 221, 
        mu_lim, r_lim, vm, nu_flatten, bar)
    fig.add_subplot(sp0)
    if data1 != None:
        sp1 = single_stripe_mur(data1, wedge, mag, scale, color, 222, 
            mu_lim, r_lim, vm, nu_flatten, bar=0)
        fig.add_subplot(sp1)
    if data2 != None:
        sp2 = single_stripe_mur(data2, wedge, mag, scale, color, 223,
            mu_lim, r_lim, vm, nu_flatten, bar=0)
        fig.add_subplot(sp2)
    if data3 != None:
        sp3 = single_stripe_mur(data3, wedge, mag, scale, color, 224, 
            mu_lim, r_lim, vm, nu_flatten, bar=0) 
        fig.add_subplot(sp3)
    if outname == None:  plt.show()
    else:  plt.savefig((outname+".ps"), papertype='letter')
    plt.close('all')

def single_stripe_mur(data, wedge, mag=0, scale=0, color=1, position=111,
                    mu_lim=None, r_lim=None, vm=None, nu_flatten=0, bar=1):
    """ change 'scale' to a string tag:  None, sqrt, flatten, log? """
    #l, b = data[:,0], data[:,1]
    x_size, y_size = 0.5, 0.5
    if mag==1:  r = coor.getr(data[:,2])
    else:       r = data[:,2]
    print len(r)
    if r_lim == None:  r_lim = (np.ma.min(r), np.ma.max(r))
    ra, dec = coor.lbToEq(data[:,0],data[:,1])
    mu, nu = coor.EqToGC(ra, dec, wedge)
    if mu_lim==None:  mu_lim = (np.ma.min(mu), np.ma.max(mu))
    x, y = r*sc.cos(mu*rad), r*sc.sin(mu*rad)
    x_bins = int((np.ma.max(x)-np.ma.min(x))/x_size) + 1
    y_bins = int((np.ma.max(y)-np.ma.min(y))/y_size) + 1
    H, x, y = np.histogram2d(y, x, [y_bins, x_bins]) #, range=[[-155.0, 110.0],[-190.0, 0.0]])
    extent = [np.ma.min(y), np.ma.max(y), np.ma.max(x), np.ma.min(x)]
    if nu_flatten == 1:
        print "!!! DOUBLE-CHECK NU FLATTENING! {0}".format(H.shape)
        for i in range(H.shape[0]):
            for j in range(H.shape[1]):
                dist = ma.sqrt( ((i-(H.shape[0]/2.0))*y_size)**2 + (j*x_size)**2)
                H[i,j] = H[i,j] / (dist*2.5)
        H = sc.sqrt(H) #sc.log10(H)
        vm = 1.5
    elif scale==1:
        H=sc.sqrt(H)
        if vm==None:  vm=16.0
    else:
        if vm==None:  vm=300.0
    """ Begin Figure """
    #fig = plt.figure(frameon=False)
    sp = plt.subplot(position)
    if color==1:  cmap = spectral_wb
    else:  cmap = 'gist_yarg'
    plt.imshow(H, extent=extent, interpolation='nearest', vmin=0.0, vmax=vm,
               cmap=cmap)
    if bar == 1:
        bar = plt.colorbar(orientation='vertical')
        if scale == 1:
            bar_ticks = sc.arange(0.0, vm+1.0, 2.0)
            bar_labels = []
            for i in range(len(bar_ticks)):   
                bar_labels.append(str(int(bar_ticks[i]**2)))
            bar.set_ticks(bar_ticks)
            bar.set_ticklabels(bar_labels, update_ticks=True)
    """ Draw axes - mu """
    mu_axis = gen_mu_axis(mu_lim, r_lim[1]) 
    plt.plot(mu_axis[0], mu_axis[1], 'k-')  # Draws main axis
    for i in range(len(mu_axis[2])):  
        x0 = [1.0*r_lim[1]*sc.cos(mu_axis[2][i]*rad), 1.1*r_lim[1]*sc.cos(mu_axis[2][i]*rad)]
        y0 = [1.0*r_lim[1]*sc.sin(mu_axis[2][i]*rad), 1.1*r_lim[1]*sc.sin(mu_axis[2][i]*rad)]
        plt.plot(x0,y0, 'k-')  # Draws mu ticks
        if mu_axis[3][i] <= 180.0:  rot = -90.0 + mu_axis[3][i]
        else:  rot =  -360.0 + mu_axis[3][i] - 90.0
        if wedge < 40:  offset = 1.12  #This accounts for different whitespace on each side of 180
        else:  offset = 1.14
        x00, y00 = offset*r_lim[1]*sc.cos(mu_axis[2][i]*rad), offset*r_lim[1]*sc.sin(mu_axis[2][i]*rad)
        #Format tick labels - this part adds whitespace to center labels on ticks
        pre, mulbl = "$", str(mu_axis[3][i])  
        if len(mulbl) < 3:  
            for k in range(3-len(mulbl)):  pre = " "+pre
        if i == (len(mu_axis[2]))-1:
            plt.text(x00, y00, "$\mu="+str(mu_axis[3][i])+"^{\circ}$", rotation=rot, 
                fontsize=12, horizontalalignment='center') #float(mu_axis[3][i]) )
        else:  plt.text(x00, y00, pre+str(mu_axis[3][i])+"^{\circ}$", rotation=rot, 
                fontsize=12, horizontalalignment='center') #float(mu_axis[3][i]) )
    """ Draw axes - r """
    r_axis = gen_r_axis(mu_lim, r_lim)
    plt.plot(r_axis[0], r_axis[1], 'k-')
    mu1 = sc.arange(mu_lim[0]-6.0, mu_lim[1]+7.0, 1.0)
    for tick in r_axis[2]:
        x1, y1 = tick*sc.cos(mu1*rad), tick*sc.sin(mu1*rad)
        if mu1[-1] < 90.0 or mu1[-1] > 360.0:
            x2, y2 = tick*sc.cos((mu1[-1])*rad), tick*sc.sin((mu1[-1])*rad)
            hl = 'right'
        else:  x2, y2 = tick*sc.cos((mu1[1])*rad), tick*sc.sin((mu1[1])*rad); hl='left'
        plt.plot(x1[:7], y1[:7], "k-")
        plt.plot(x1[-7:], y1[-7:], "k-")
        if tick == r_axis[2][-1]:
            plt.text(x2, y2, "$r="+str(int(tick))+"$", rotation=0.0, fontsize=12,
                 horizontalalignment=hl, verticalalignment='bottom')
        else:  plt.text(x2, y2, "$"+str(int(tick))+"$", rotation=0.0, fontsize=12,
                 horizontalalignment=hl, verticalalignment='bottom')
    # Draw axes - g (TO BE DONE)
    """ Clean up output """
    y_lim = plt.ylim()
    plt.ylim(y_lim[1], y_lim[0])
    plt.setp(sp.get_xticklabels(), visible=False)
    plt.setp(sp.get_yticklabels(), visible=False)
    plt.setp(sp.get_xticklines(),  visible=False)
    plt.setp(sp.get_yticklines(),  visible=False)
    #ax = fig.add_axes([0, 0, 1, 1])  #These two lines supress the axes... supposedly
    #ax.axis('off')
    return sp

def gen_mu_axis(mu_lim, edge):
    d = edge*1.05
    low, high = ma.floor(mu_lim[0]/15.0), ma.ceil(mu_lim[1]/15.0)
    #low, high = (mu_lim[0]/15.0), (mu_lim[1]/15.0)
    mu = sc.arange((low*15.0)-10.0, (high*15.0)+10.1, 1.0)
    x, y = d*sc.cos(mu*rad), d*sc.sin(mu*rad)
    mu_ticks = sc.arange(low, high+1.0, 1.0)
    mu_ticks = mu_ticks*15.0
    #u, v = d*sc.cos(mu_ticks*rad), d*sc.sin(mu_ticks*rad)
    mu_labels = []
    for tick in mu_ticks:
        if tick >= 360.0:
            mu_labels.append(int((tick-360.0) ) )
        else:  mu_labels.append(int(tick) )
    return [x,y,mu_ticks,mu_labels]
    
def gen_r_axis(mu_lim, r_lim):
    low, high = mu_lim[0]-3.0, mu_lim[1]+3.0
    x, y = [], []
    r = sc.array([1.05*r_lim[1], 0.0, 1.05*r_lim[1]])
    mu = sc.array([low, 0.0, high])
    x, y = r*sc.cos(mu*rad), r*sc.sin(mu*rad)
    r_ticks = sc.arange(10.0, 10.0*ma.floor(r_lim[1]/10.0)+1.0, 10.0)
    return [x,y,r_ticks]

def plot_stripe_mu(data, wedge, musize=1.0, name="out"):  #Aspect ratio is crap - just bin along mu?
    l, b, r = data[:,0], data[:,1], coor.getr(data[:,2])
    ra, dec = coor.lbToEq(l,b)
    mu, nu = coor.EqToGC(ra, dec, wedge)
    mu_hist = np.histogram(mu, int((np.ma.max(mu)-np.ma.min(mu))/musize), (np.ma.min(mu),np.ma.max(mu)))
    plt.figure(1)
    plt.bar(mu_hist[1][:-1], mu_hist[0], musize)
    #plt.savefig((name+".ps"), papertype='letter')
    plt.show()
    plt.close('all')

def plot_wedge_density(data, wedge=82, q=0.458, r0=19.5, mu_min=310.0, mu_max=419.0,
                       streams=None, perturb=None, name="out", mag=0, plot=1, size=1.0):
    """Plots the average density of a wedge versus radius
    stream coords. should be a list of streams, with each sub-list being a
    standard 6-parameter set of stream parameters; perturb is weight of
    perturbation, if any;  mag is whether data is in g magnitudes or not.
    CURRENTLY ONLY WORKS (WELL) FOR STRIPE 82!!!"""
    #Get Average Density with r
    x, y, z = coor.lbr2xyz(data[:,0], data[:,1], data[:,2])
    if mag==1: r = coor.getr(data[:,2])
    else:  r = data[:,2]
    max, min = np.ma.max(r), np.ma.min(r)
    bins = int((max - min + 1.0) / size)
    rho_raw = sc.zeros(bins, int)
    for i in range(len(r)):
        index = int( (r[i]-min) / size)
        rho_raw[index] = rho_raw[index] + 1
    # get volume at each distance
    #mu_min, mu_max = 310.0, 419.0  #Defined in def inputs
    nu_min, nu_max = (-1.25*ma.pi/180.0), (1.25*ma.pi/180.0)
    mu_bit = mu_max - mu_min
    nu_bit = sc.sin(nu_max) - sc.sin(nu_min)
    V = sc.zeros(bins, float)
    for i in range(bins):
        r_bit = (min + (i+1)*size)**3 - (min + i*size)**3
        V[i] = r_bit*(1./3.)*mu_bit*nu_bit
    # make the histogram of rho(r)
    rho_hist = sc.zeros((bins,3), float)
    for i in range(bins):
        rho_hist[i,0] = min + (i*size)
        rho_hist[i,1] = rho_raw[i] / V[i]
        rho_hist[i,2] = (ma.sqrt(rho_raw[i] + 1.0) + 1.0) / V[i]
    """Plot Herquist, stream, perturbation functions over histogram"""
    #get gal-centered herquist distribution for r points in wedge
    mu_avg = ((mu_max - mu_min) / 2.0) + mu_min
    nu_avg = 0.0
    sol_r = sc.arange(min, max, size)
    x,y,z = coor.GC2xyz(mu_avg, nu_avg, sol_r, wedge)
    Hern_rho = twg.hernquist_profile(x,y,z,q,r0)
    H_scale = sc.sum(rho_hist[:,1]) / sc.sum(Hern_rho) #multiply a factor of 2 to account for different # of bins
    for i in range(len(Hern_rho)): Hern_rho[i] = Hern_rho[i]*H_scale
    # DO STREAM!!!
    ################
    #get gal-centered perturbation distribution for r points in wedge
    if perturb != None:
        perturb_rho = twg.perturb_density(x,y,z, (0.0, 0.79, -19.9) )
        p_scale = 1.0
        for i in range(len(perturb_rho)): perturb_rho[i] = perturb_rho[i]*p_scale
    #Plot it all
    if plot == 1:
        fig = plt.figure()
        plt.bar(rho_hist[:,0], (rho_hist[:,1]/Hern_rho)-1.0, size)
        #plt.plot(sol_r, Hern_rho, c='r')
        if perturb != None:
            plt.plot(sol_r, perturb_rho, c='r')
        plt.xlabel(r'Sun-centered $r$ (kpc)')
        plt.ylabel(r'$\rho$ (stars/kpc)')
        #plt.savefig((name+".ps"), papertype='letter')
        plt.show()
        plt.close('all')
    return rho_hist

def plot_subwedge_density():
    return 0

def plot_double_mur(data, wedge=82, mu_min=310.0, mu_max=419.0, name="out",
                    mag=0, rsize=2.0, musize=2.0):
    """ !!! TAKE INTO ACOUNT mu WRAPAROUND AT 360.0 !!!"""
    ra, dec = coor.lbToEq(data[:,0],data[:,1])
    mu, nu = coor.EqToGC(ra, dec, wedge)
    r_min, r_max = int(sc.ma.min(data[:,2])), int(sc.ma.max(data[:,2]))+1
    r_hist = np.histogram(data[:,2], int((r_max-r_min)/rsize), (r_min, r_max))
    mu_hist = np.histogram(mu, int((mu_max-mu_min)/musize), (mu_min, mu_max))
    # plot r hist
    fig = plt.figure(1)
    plt.subplot(211)
    plt.bar(r_hist[1][:-1], r_hist[0], rsize)
    plt.xlabel("r, kpc")
    # plot mu hist
    plt.subplot(212)
    plt.bar(mu_hist[1][:-1], mu_hist[0], musize)
    plt.xlabel("mu")
    plt.show()
    plt.close('all')

def do_compare_wedges(file1="stars-82.txt", file2="Stripe82_coadd.csv", stripe=82, hern=0):
    """ Modify if size is not 1.0 """
    one_run = fi.read_data(file1)
    or_l = len(one_run[:,0])
    or_hist = plot_wedge_density(one_run, stripe, q=0.458, r0=19.4, streams=None, perturb=None,
                       name="_rho1", mag=0, plot=0)
    coadd = fi.read_data(file2)
    ca_l = len(coadd[:,0])
    ca_hist = plot_wedge_density(coadd, stripe, q=0.458, r0=19.4, streams=None, perturb=None,
                       name="_rho2", mag=0, plot=0)
    # Separate into heights
    #print or_hist[:,0]
    #print ca_hist[:,0]
    or_h = or_hist[:,1]
    ca_h = ca_hist[:,1]
    if hern==1:
        # Get Hernquist profile
        mu_avg = ((419.0 - 310.0) / 2.0) + 310.0
        nu_avg = 0.0
        sol_r = sc.arange( (np.ma.min(ca_hist[:,0])+0.5), (np.ma.max(ca_hist[:,0])+1.5), (1.0))
        x,y,z = coor.GC2xyz(mu_avg, nu_avg, sol_r, 82)
        Hern_rho = twg.hernquist_profile(x,y,z,q=0.458,r0=19.4)
        H_scale = sc.ma.max(Hern_rho)
        for i in range(len(Hern_rho)): Hern_rho[i] = Hern_rho[i]/H_scale
        for i in range(len(or_h)):
            or_h[i] = or_h[i]/Hern_rho[i]
        norm = sc.ma.max(ca_h)
        for i in range(len(ca_h)):
            ca_h[i] = ca_h[i]/Hern_rho[i]
        # Differentiate the bins
        if len(or_h) > len(ca_h):  l = len(or_h)
        else:  l = len(ca_h)
        diff_h = sc.zeros(l)
        for i in range(len(or_h)):
            diff_h[i] = diff_h[i] + or_h[i]
        for i in range(len(ca_h)):
            diff_h[i] = diff_h[i] - ca_h[i]
    else:
        # Divide the first data set by the second
        if len(or_h) < len(ca_h):
            l = len(or_h)
            extra_h = -0.1*sc.ones((len(ca_h)-l))
        else:
            l = len(ca_h)
            extra_h = 0.1*sc.ones((len(or_h)-l))
        diff_h = sc.zeros(l)
        for i in range(l):
            #diff_h[i] = ( (or_h[i]/or_l) / (ca_h[i]/ca_l) ) - 1.0
            diff_h[i] = ( or_h[i] / ca_h[i] ) - 1.0
    #Make plot
    fig = plt.figure()
    #plt.bar(ca_hist[:,0], ca_hist[:,1], 1.0, fill=None, ec="b")
    #plt.bar(or_hist[:,0], or_hist[:,1], 1.0, fill=None, ec="r")
    #plt.ylabel(r'$\rho$ (stars/kpc)')
    plt.bar(ca_hist[:l,0], diff_h, 1.0)
    if hern != 1:
        #plt.bar(ca_hist[l:,0], extra_h, 1.0, ec="r", fill=None)
        plt.plot([np.ma.min(ca_hist[:,0]),np.ma.max(ca_hist[:,0])], [0.0, 0.0], "k:")
    plt.xlabel(r'Sun-centered $r$ (kpc)')
    plt.ylabel(r'$(\rho_{unmatched} / \rho_{matched}) - 1.0$ (stars/kpc)')
    #plt.savefig(("diff.ps"), papertype='letter')
    plt.show()
    plt.close('all')

def do_wedgeplots():
    files = glob.glob('stripe_[7-8][0-9].txt')
    print files
    for stripe in files:
        wedge = stripe[7:9]
        data = fi.read_data(stripe)
        plot_stripe_mur(data, int(wedge), ("wedge_"+wedge))

def do_testwedgeplots():
    files = glob.glob('single*.txt')
    print files
    for stripe in files:
        wedge = stripe[7:-4]
        data = fi.read_data(stripe)
        #plot_stripe_mur(data, 82, ("wedge_"+wedge))
        plot_stripe_mug(data, 82, ("wedge_g_"+wedge))

def do_densityplots():
    files = glob.glob('single*[0-4]_2.txt')
    #files = glob.glob('stripe_[7-8][0-9].txt')
    print files
    for stripe in files:
        #wedge = stripe[7:9]
        wedge = stripe[7:12]
        data = fi.read_data(stripe)
        plot_wedge_density(data, wedge=82, q=1.0, r0=19.5, mu_min=310.0, mu_max=419.0,
                       streams=None, perturb=None, name=("hernquist_"+wedge), mag=0) #int(wedge)

def do_onewedgeplot(filename="no_matches_82.txt", wedge=82, delimiter=None):
    data = fi.read_data(filename, delimiter)
    plot_stripe_mur(data, wedge, filename[:-4], mag=0)

    
if __name__ == "__main__":
    if len(sys.argv) > 1: data = fi.read_data(sys.argv[1])
    if len(sys.argv) > 2: wedge = int(sys.argv[2])
    else:  wedge = 82
    #plot_stripe_lb(data)
    #plot_stripe_3D(data)
    #plot_stripe_mur(data, wedge, mag=0, scale=1) #, mu_lim=(310.0, 419.0))
    #plot_stripe_mug(data, wedge)
    #plot_stripe_mu(data, wedge)
    #plot_stripe_munu(data, wedge)
    #plot_wedge_density(data, wedge, q=0.458, r0=19.4, name=sys.argv[1][:-4]+"_rho", mag=0)
    #plot_double_mur(data, wedge=82, mu_min=310.0, mu_max=419.0,
    #               name=sys.argv[1][:-4]+"_mur", mag=0, rsize=1.0, musize=1.0)
    #do_wedgeplots()
    #do_testwedgeplots()
    #do_densityplots()
    #do_onewedgeplot("stripe-10-DR7-clean.txt", 10, None)
    #do_compare_wedges(file1="unmatched-82-1arcsec.txt", file2="matched-82-1arcsec.txt", stripe=82, hern=0)
    lbpolar_plot("./sep_lbr/", hemi='N', bin_size=0.5, color=1, scale=1, primary=1)
    # /home/newbym2/Desktop/star_holder/sep_lbr
    
    """t0 = time.time()
    wedges = range(70,90)
    all_stripes(wedges, "FTO_south_mnewby.csv")
    tf = time.time()
    print "Total time elapsed: {0} minutes".format(((tf-t0)/60.0))"""
    
    print "# --- Done"
