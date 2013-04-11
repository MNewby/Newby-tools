import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import gradient_descent as gd
import functions as func
import poly_fit as poly
import monte_carlo as mc

'''python script for fitting Girardi isochrones to CMDs, and getting errors.
Matthew Newby, September 18, 2011'''

def format_isochrone(iso_in, cuts=(8.0, 3.5)):
    Mg_raw, gminusr_raw = iso_in[:,8], (iso_in[:,8]-iso_in[:,9])
    for i in range(len(Mg_raw)):
        if Mg_raw[i] < cuts[0]:  a = i; break
    for i in range(len(Mg_raw)):    
        if Mg_raw[i] < cuts[1]:  b = i; break
    Mg, gminusr = Mg_raw[a:b], gminusr_raw[a:b]
    for i in range(len(Mg)):  #correction
        gminusr[i] = gminusr[i] -0.015*Mg[i] + 0.089
    return sc.array(zip(Mg,gminusr, 0.01*sc.ones(len(Mg))))
    
def format_data(data_in, dist, cuts=(3.5, 8.0), bin_size=0.2):  #cuts=(16.0, 24.5)
    """ Takes in raw data, bins it, then finds averages and sigmas through
        2.5-sigma fit???"""
    g0 = data_in[:,2] - 5.0*(sc.log10(dist*1000.0) - 1.0)
    gminusr = (data_in[:,2]-data_in[:,3])
    # Get bin edges
    edges = sc.arange(cuts[0], (cuts[1]+bin_size), bin_size)
    centers = edges[:-1] + (bin_size/2.0)
    # Get bin averages, stdevs
    middles, averages, stdevs = [], [], []
    for i in range(len(edges)-1):
        low, high = edges[i], edges[i+1]
        gr_holder = []
        for j in range(len(g0)):
            if (g0[j] > low) and (g0[j] <= high):  gr_holder.append(gminusr[j])
        avg, std = rejection_fit(gr_holder)
        if  std == 0.0:   continue
        middles.append(centers[i])
        averages.append(avg)
        stdevs.append(std)
    return sc.array(zip(centers, averages, stdevs))

def rejection_fit(data_in, threshold=2.5):
    holder, correct_exit = [], 0
    for i in range(len(data_in)):  holder.append(data_in[i])  #make it local
    while len(holder) > 1:
        avg = sc.mean(holder)
        std = sc.std(holder)
        out = []
        for i in range(len(holder)):
            if abs(avg - holder[i]) <= threshold*std:
                out.append(holder[i])
        if len(holder) == len(out):  break
        holder = []
        for i in range(len(out)):  holder.append(out[i])
    if len(holder) < 2:
        return (None, 0.0)
    else:
        return sc.mean(out), sc.std(out)

def R_squared(clus, params):
    RR = 0
    y = func.Nth_order_polynomial(clus[:,0], params)
    for i in range(len(clus[:,0])):
        top = (clus[i,1]-y[i])*(clus[i,1]-y[i])
        bottom = clus[i,2]*clus[i,2]
        RR = (top/bottom) + RR
    return RR

def fit_iso(name, dist, plot=1):
    if name != 'Pal5':
        clus_data = fi.read_data("noU_NGC_"+name+"_cluster.csv", ",")
    else:
        clus_data = fi.read_data("noU_"+name+"_cluster.csv", ",")
    if dist < 10.0:  high = 7.0
    elif dist < 20.0:  high = 6.0
    else:  high = 5.5
    clus = format_data(clus_data, dist, cuts=(3.5, high))
    iso_data, age = [], []
    for i in range(len(mod)):
        try:
            filename = "./Giso_shifted/"+name+mod[i]+".dat"
            iso_data.append(fi.read_data(filename))
            age.append(d_age[i])
        except IOError:
            print "!!! File not found - {0}".format(filename)
            continue
    """ Get R-squared values for each isochrone"""
    iso, RR = [], []
    for i in range(len(iso_data)):
        iso.append(format_isochrone(iso_data[i], cuts=(7.0, 3.5)))
        results = poly.poly_fit(0.0,iso[i][:,0],iso[i][:,1],iso[i][:,2], 6, verbose=0)
        RR.append(R_squared(clus, results))
    points = sc.array(zip(age, RR))
    max = sc.ma.max(points[:,1])
    RRmod = -1.0*(points[:,1] - max)
    sigma = 0.1*sc.ones(len(RRmod))
    c, d, start = mc.do_MCMC(func.gaussian_function, sc.array([0.1, 0.0, 0.2]), sc.array([0.001,0.001,0.001]),
               points[:,0], RRmod, sigma, "test", number_steps=10000, save=0)
    best = gd.gradient_descent(func.gaussian_function, start,
                               points[:,0], RRmod, sigma)
    # One-dimensional errors
    if len(RR) != 5:
        error = np.sqrt(abs((8.0*0.2*0.2) / (2.0*(RR[-1]-RR[0])) ) )  #Uses last point, which should be the '-2'
    else:
        error = np.sqrt(abs((8.0*0.2*0.2) / (RR[2] + RR[4] - 2.0*RR[0])))  # Hesssian error for one parameter fit
    # Plot these bad boys
    if plot==1:
        plt.figure(1)
        plt.subplot(211)
        plt.errorbar(clus[:,0], clus[:,1], yerr=clus[:,2], fmt='o')
        for i in range(len(iso)):
            plt.plot(iso[i][:,0], (iso[i][:,1]) )
        plt.subplot(212)
        plt.scatter(points[:,0], points[:,1])
        x = sc.arange(plt.xlim()[0], plt.xlim()[1], 0.01)
        y = -1.0*(func.gaussian_function(x, best)) + max #func.gaussian_function(x, best) #
        plt.plot(x,y, 'g-')
        plt.savefig(name+"_isoAn2.ps", papertype='letter')
        plt.close('all')
    return (sc.array(zip(age,RR)), best, error)

def fit_distance():
    dmod = [0.0, 5.0, 10.0, -5.0, -10.0]
    filename = "./Girardi_Isochrones/"+name+"_0.dat"
    iso_data = fi.read_data(filename)
    if name != 'Pal5':
        clus_data = fi.read_data("noU_NGC_"+name+"_cluster.csv", ",")
    else:
        clus_data = fi.read_data("noU_"+name+"_cluster.csv", ",")
    if dist < 10.0:  high = 7.0
    elif dist < 20.0:  high = 6.0
    else:  high = 5.5
    clus = format_data(clus_data, dist, cuts=(3.5, high))


if __name__ == "__main__":
    name = ['4147', '5024', '5053', '5272', '5466', '5904', '6205', '6341', '7078', '7089', 'Pal5']
    distance = [19.3, 18.7, 18.5, 10.4, 15.6, 8.0, 7.7, 8.7, 11.0, 11.5, 21.0]
    #modulus = 5.0*(sc.log10(distance*1000.0) - 1.0)
    skip = [0,0,0,0,0,0,0,0,0,0,0]  #set to 1 to skip a cluster
    mod = ['_0', '_1', '_2', '-1', '-2']
    d_age = [0.0, 0.2, 0.4, -0.2, -0.4]
    # Fit these guys
    stuff, err = [], []
    for i in range(len(name)):
        if skip[i] == 1:  continue
        print "###-- Starting run: {0} from {1} kpc".format(name[i], distance[i])
        goods, params, error = fit_iso(name[i], distance[i], plot=0)
        #print "$$$", goods, params, error
        stuff.append(params)
        err.append(error)
    print "Amplitude - Mean Offset - Sigma - Hessian error"
    for i in range(len(stuff)):
        #print "{0} {1} {2} {3}".format(stuff[i][0], stuff[i][1], stuff[i][2], err[i])
        print err[i]
    print "#--Done"
    

    """
    clus_data = fi.read_data("noU_NGC_6205_cluster.csv", ",")
    clus = format_data(clus_data, 8.0, cuts=(3.5, 7.0))
    iso_data = fi.read_data("./Girardi_Isochrones/6205_0.dat")
    iso = format_isochrone(iso_data, cuts=(8.0, 3.5))
    results = poly.poly_fit(0.0,iso[:,0],iso[:,1],iso[:,2], 6)
    x = sc.arange(3.5, 8.0, 0.1)
    y = func.Nth_order_polynomial(x, results)
    #fi.write_data(iso, "test.txt")
    plt.figure(1)
    plt.plot(iso[:,0], (iso[:,1]) )
    plt.errorbar(clus[:,0], clus[:,1], yerr=clus[:,2], fmt='o')
    plt.plot(x,y, 'r:')
    plt.show()
    """
"""

d0 = fi.read_data("./Girardi_Isochrones/5024_0.dat")
d1 = fi.read_data("./Girardi_Isochrones/5024_1.dat")
d2 = fi.read_data("./Girardi_Isochrones/5024_2.dat")
d3 = fi.read_data("./Girardi_Isochrones/5024-1.dat")
d4 = fi.read_data("./Girardi_Isochrones/5024-2.dat")

plt.figure(1)
plt.plot(d0[:,8]-d0[:,9], d0[:,8])
plt.plot(d1[:,8]-d1[:,9], d1[:,8])
plt.plot(d2[:,8]-d2[:,9], d2[:,8])
plt.plot(d3[:,8]-d3[:,9], d3[:,8])
plt.plot(d4[:,8]-d4[:,9], d4[:,8])
plt.xlim(0.0,0.4)
plt.ylim(10.0, 0.0)
plt.show()


names = ['4147', '5024', '5053', '5272', '5466', '5904', '6205', '6341',
         '7078', '7089', 'Pal5']

"""