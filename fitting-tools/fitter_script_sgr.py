import sys
sys.path.insert(0, '../utilities')
import glob
import math as ma
import numpy as np
import scipy as sc
import functions as func
import files as fi
import fit as fit
import matplotlib
import matplotlib.pyplot as plt

path = "/home/newbym2/Dropbox/Research/sgrLetter/hist_txts"
files = glob.glob(path+"/*.out")
print files

for f in files:
    # load in data
    name = f[-9:-4]
    data = np.loadtxt(f)
    for i in range(len(data[:,0])):
        data[i,1] = data[i,1] - 200.0
        if data[i,1] < 0.0:  data[i,1] = 0.0
    # Set up data
    x, y = data[:,0], data[:,1]
    e = func.poisson_errors(y)
    #fit it
    fitter = fit.ToFit(x,y,e)
    #fitter.function=func.double_gaussian
    #fitter.update_params([20.0, 10.0, 5.0, 20.0, 0.0, 10.0])
    #fitter.step = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    #fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma"]
    fitter.function=func.quad_fat_gaussians
    fitter.update_params([20.0, 10.0, 5.0, 20.0, 0.0, 5.0, 5.0, 15.0, 5.0, 15.0])
    fitter.step = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma", "amp","sigma", "amp","sigma"]
    print "\n {0}:".format(name)
    print "# - Starting from:", fitter.params, fitter.RR
    print "# - With steps:", fitter.step
    path1 = fit.gradient_descent(fitter, its=10000, line_search=0)
    path2 = fit.MCMC(fitter)
    #errors = fit.get_Hessian_errors(fitter)
    # cut from ToFit
    new_params=fitter.params
    plt.figure(1)
    plt.bar(fitter.x, fitter.y, ((np.ma.max(fitter.x)-np.ma.min(fitter.x))/len(fitter.x)), align='center')
    plt.plot(fitter.x, fitter.function(fitter.x, new_params), 'k-')
    print new_params[:4]
    plt.plot(fitter.x, func.gaussian_function(fitter.x, new_params[:3]), 'k--')
    plt.plot(fitter.x, func.gaussian_function(fitter.x, new_params[3:]), 'k--')
    plt.xlabel("B")
    plt.ylabel("counts")
    plt.title(name, fontsize=8 )
    out_name = name+"_plot_"+".png"
    #plt.show()
    plt.savefig(out_name)
    plt.close('all')
    # counter allows for a series of uniquely-named plots


#fitter.range = [(0, 0.0, 5.1, 0.25), (1, 80.0, 181.0, 1.0), (2, 5.0, 55.5, 0.5), (3, 0.0, 5.1, 0.05)]
#yi, xi = 1, 2
#surface = fit.get_surface(fitter, yi, xi, plot=1, path=(path1[:,yi], path1[:,xi]))
#surface = fit.get_surface(fitter, yi, xi, plot=1, path=(path2[:,yi], path2[:,xi]))
#fitter.plot_type = "bar"
#fitter.plot()

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
