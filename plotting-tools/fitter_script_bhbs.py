import math as ma
import numpy as np
import scipy as sc
import functions as func
import files as fi
import fit as fit
import matplotlib
import matplotlib.pyplot as plt

run = range(1,6)
start_params = [3.0, 17.0, 0.75, 5.0]
bins=20
#run=[1]  #debug

for i in run:
    print "# --- Fitting region {0}".format(i)
    data = fi.read_data("../Cetus_Polar/bhb_dmod"+str(i)+".txt", skip=1)
    print "# --- l-coordinate mean, std: {0}, {1}".format(np.mean(data[:,0]), np.std(data[:,0]) )
    b_hist, b_edges = np.histogram(data[:,1], bins, range=[80.0, 180.0])
    b_err = func.poisson_errors(b_hist)
    b_wid = (np.ma.max(data[:,1]) - np.ma.min(data[:,1]) ) / bins
    d_hist, d_edges = np.histogram(data[:,2], bins, range=[16.0, 19.0])
    d_err = func.poisson_errors(d_hist)
    d_wid = (np.ma.max(data[:,2]) - np.ma.min(data[:,2]) ) / bins
    # Initialize
    fitter = fit.ToFit(d_edges[:-1]+(d_wid/2.0), d_hist , d_err)
    fitter.function=func.gauss_plus_floor
    fitter.update_params(start_params)
    fitter.step = [0.01, 0.01, 0.01, 0.01]
    fitter.param_names = ["amp", "mu", "sigma", "floor"]
    fitter.min = [None, None, 0.0, None]
    fitter.max = [None, None, None, None]
    print "# - Starting from:", fitter.params, fitter.RR
    print "# - With steps:", fitter.step
    # Fit
    path1 = fit.gradient_descent(fitter, its=10000, line_search=0)
    errors = fit.get_Hessian_errors(fitter)
    gd_params = fitter.params, fitter.RR
    #fitter.update_params(start_params)
    print "# - Starting from:", fitter.params, fitter.RR
    path2 = fit.MCMC(fitter, iters=10000, annealing=0, verbose=0)
    mc_params = fitter.params, fitter.RR
    # Surfaces
    fitter.range = [(0, 3.5, 4.5, 0.02), 
                    (1, 17.0, 18.0, 0.02), 
                    (2, 0.0, 0.4, 0.008), 
                    (3, 0.0, 10.0, 0.1)]
    yi, xi = 1, 2
    surface = fit.get_surface(fitter, yi, xi, plot="surface"+str(i)+"_big", path=(path2[:,yi], path2[:,xi]))
    # Plot of statistics
    plt.figure()
    plt.subplot(211)
    plt.hist(path2[:,1], 50, color="grey")
    plt.hist(path2[5000:,1], 50, color="white")
    plt.xlabel(r"$\mu$")
    plt.subplot(212)
    plt.hist(path2[:,2], 50, color="grey")
    plt.hist(path2[5000:,2], 50, color="white")
    plt.xlabel(r"$\sigma$")
    plt.savefig("region"+str(i)+"_dist_big.png", papertype="letter")
    #plt.show()
    plt.close('all')
    #Plots 
    plt.figure(i)
    plt.bar(d_edges[:-1], d_hist, d_wid, color="grey")
    if gd_params[1] < mc_params[1]:  best = gd_params[0]
    else:  best = mc_params[0]
    x = np.arange(16.0, 19.0, 0.01) 
    y = func.gauss_plus_floor(x, best)
    plt.plot(x,y,"k-")
    #plt.xlim(80.0, 180.0)
    plt.savefig("region"+str(i)+"_big.png", papertype="letter")
    #plt.show()
    plt.close("all")
    
    
    
    

"""  ------------------------  OLD ----------------------- 
run = 5
print "# --- Fitting region b{0}".format(run)

data = fi.read_data("bhb_histogram.txt", ",")
back = fi.read_data("bhb_histogram_wide.txt", ",")

x = data[:,0]
y = 2.0*data[:,run] - back[:,run]

x = back[:,0]

s = back[:,run+5]
e = func.poisson_errors(back[:,run]/s)*s

fitter = fit.ToFit(x,y,e)
fitter.function=func.gaussian_function
fitter.update_params([3.0, 136.0, 30.0]) #, 2.0])
fitter.step = [1.0, 10.0, 2.5] #, 2.0]
fitter.param_names = ["amp", "mu", "sigma", "floor"]
print "# - Starting from:", fitter.params, fitter.RR
print "# - With steps:", fitter.step

path1 = fit.gradient_descent(fitter, its=10000, line_search=0)
errors = fit.get_Hessian_errors(fitter)

fitter.update_params([3.0, 140.0, 42.0, 1.0])
path2 = fit.MCMC(fitter, iters=10000, annealing=0, verbose=0)

fitter.range = [(0, 0.0, 5.1, 0.25), (1, 80.0, 181.0, 1.0), (2, 5.0, 55.5, 0.5), (3, 0.0, 5.1, 0.05)]

yi, xi = 1, 2
surface = fit.get_surface(fitter, yi, xi, plot=1, path=(path1[:,yi], path1[:,xi]))
#surface = fit.get_surface(fitter, yi, xi, plot=1, path=(path2[:,yi], path2[:,xi]))

fitter.plot_type = "bar"
fitter.plot()

# Plot of statistics
plt.figure()
plt.subplot(211)
plt.hist(path2[:,1], 50)
plt.subplot(212)
plt.hist(path2[:,2], 50)
#plt.scatter(path2[:,0], path2[:,1])
plt.show()
"""
