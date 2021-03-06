import sys
sys.path.insert(0, '../utilities')
import math as ma
import numpy as np
import scipy as sc
import functions as func
#import files as fi
import fit as fit
import matplotlib
import matplotlib.pyplot as plt

#run = 5

#print "# --- Fitting region b{0}".format(run)

data = np.loadtxt("/home/newbym2/Dropbox/Research/sgrLetter/virgo_edge_points.txt")
x = data[:-3,0]
y = data[:-3,1]
e = func.poisson_errors(y)

#back = fi.read_data("bhb_histogram_wide.txt", ",")
#x = back[:,0]
#s = back[:,run+5]

fitter = fit.ToFit(x,y,e)
fitter.function=func.get_2gauss_y
fitter.update_params([14.0, 250.0, 30.0, 100.0])
fitter.step = [1.0, 1.0, 1.0, 1.0]
fitter.param_names = ["amp", "mu", "sigma", "floor"]
print "# - Starting from:", fitter.params, fitter.RR
print "# - With steps:", fitter.step

path1 = fit.gradient_descent(fitter, its=10000, line_search=0)
errors = fit.get_Hessian_errors(fitter)

#fitter.update_params([3.0, 140.0, 42.0, 1.0])
path2 = fit.MCMC(fitter, iters=10000, annealing=0, verbose=0)

fitter.range = [(0, 0.0, 5.1, 0.25), (1, 80.0, 181.0, 1.0), (2, 5.0, 55.5, 0.5), (3, 0.0, 5.1, 0.05)]

#yi, xi = 1, 2
#surface = fit.get_surface(fitter, yi, xi, plot=1, path=(path1[:,yi], path1[:,xi]))
#surface = fit.get_surface(fitter, yi, xi, plot=1, path=(path2[:,yi], path2[:,xi]))
#fitter.plot_type = "bar"
fitter.plot()

"""
# Plot of statistics
plt.figure()
plt.subplot(211)
plt.hist(path2[:,1], 50)
plt.subplot(212)
plt.hist(path2[:,2], 50)
#plt.scatter(path2[:,0], path2[:,1])
plt.show()
"""
