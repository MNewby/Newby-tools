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

#path = "/home/newbym2/Dropbox/Research/sgrLetter/hist_fits_spec_cut"
path = "/home/newbym2/Dropbox/Research/sgrLetter/hist_fits_smallbins"
#path = "/home/newbym2/Dropbox/Research/sgrLetter/hists_r30_cutoff"
files = glob.glob(path+"/*.out")
print files

def clean_data(data):
    mask = sc.zeros(len(data[:,0]), dtype=bool)
    for i in range(len(data[:,0])):
        #data[i,1] = data[i,1] - 30.0  #spec_cut
        #data[i,1] = data[i,1] - 100.0 # all data
        #data[i,1] = data[i,1] - (1.6*data[i,0] +140.0)
        if data[i,1] <= 0.0:  mask[i] = 1  #flag zeros
        #if data[i,0] > 15.0:  mask[i] = 1 #flag virgo-ish bits in flat bg subtraction
    # make new array that contains only mask==0
    data_out = []
    for i in range(len(data[:,0])):
        if mask[i]==0:  data_out.append(data[i,:])
    return np.array(data_out)

for f in files:
    # load in data
    #name = f[-18:-13]  #spec_cut
    name = f[-9:-4]  #all_data
    data = clean_data(np.loadtxt(f) )
    if len(data)<2:  continue
    # Set up data
    x, y = data[:,0], data[:,1]
    e = func.poisson_errors(y)
    # make a line, using the average of 10 bins on each end as the anchor points
    x0, y0, x1, y1 = [], [], [], []
    i = 0
    while len(x0) < 5:
        if y[i] > 99.0:  x0.append(x[i]);  y0.append(y[i])
        i += 1
    if eval(name) > 291:  i = -20
    elif eval(name) > 274:  i = -40
    elif eval(name) < 201:  i = -20
    else:  i = -1
    while len(x1) < 5:
        if y[i] > 99.0:  x1.append(x[i]);  y1.append(y[i])
        i -= 1
    xi, yi = np.mean(x0), np.mean(y0)
    xf, yf = np.mean(x1), np.mean(y1)
    #xi, yi = np.mean(x[2:6]), np.mean(y[2:6])
    #xf, yf = np.mean(x[-7:-3]), np.mean(y[-7:-3])
    aa = (yf-yi)/(xf-xi)
    bb = yf - (aa*xf)
    #fit it
    fitter = fit.ToFit(x,y,e)
    """
    fitter.function=func.double_gauss_line
    fitter.update_params([6.0, 0.0, 5.0, 6.0, -10.0, 5.0, aa, bb])
    fitter.step = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]
    fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma", "slope", "intercept"]
    """
    fitter.function=func.quad_fat_gauss_line
    fitter.update_params([6.0, 5.0, 1.0, 6.0, -5.0, 1.0, 5.0, 10.0, 5.0, 10.0, aa, bb])
    fitter.step = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]
    fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma", "amp","sigma", "amp","sigma","slope","intercept"]
    """
    fitter.function=func.double_gaussian
    fitter.update_params([20.0, 0.0, 5.0, 20.0, 20.0, 15.0])
    fitter.step = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma"]

    fitter.function=func.quad_fat_gaussians
    fitter.update_params([6.0, 10.0, 1.0, 6.0, -10.0, 1.0, 5.0, 10.0, 5.0, 10.0])
    fitter.step = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    fitter.param_names = ["amp", "mu", "sigma", "amp", "mu", "sigma", "amp","sigma", "amp","sigma"]
    """
    print "\n {0}:".format(name)
    print "# - Starting from:", fitter.params, fitter.RR
    print "# - With steps:", fitter.step
    path1 = fit.gradient_descent(fitter, its=10000, line_search=0)
    path2 = fit.MCMC(fitter)
    #errors = fit.get_Hessian_errors(fitter)
    # cut from ToFit
    new_params=fitter.params
    plt.figure(1)
    plt.bar(fitter.x, fitter.y, ((np.ma.max(fitter.x)-np.ma.min(fitter.x))/len(fitter.x)), align='center', color='w')
    plt.plot(fitter.x, fitter.function(fitter.x, new_params), 'k-')
    print new_params[:4]
    plt.plot(fitter.x, func.gaussian_function(fitter.x, new_params[:3]), 'k--')
    plt.plot(fitter.x, func.gaussian_function(fitter.x, new_params[3:6]), 'k--')
    plt.plot(fitter.x, func.gaussian_function(fitter.x, [new_params[6], new_params[1], new_params[7]]), 'k--')
    plt.plot(fitter.x, func.gaussian_function(fitter.x, [new_params[8], new_params[4], new_params[9]]), 'k--')
    plt.xlabel("B")
    plt.ylabel("counts")
    plt.title(name, fontsize=8 )
    out_name = name+"_plot"+".png"
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
