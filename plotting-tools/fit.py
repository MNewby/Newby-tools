import math as ma
import numpy as np
import scipy as sc
import functions as func
import files as fi
import time as time
import matplotlib
import matplotlib.pyplot as plt
import Hessian_errors as he

class ToFit:
    """ A class that contains the relevant info for fitting a function to data """
    def __init__(self, data_x, data_y, data_sig):
        """ Default Parameters """
        # important values
        self.x = data_x
        self.y = data_y
        # "p" uses Poisson Errors, A stringified float uses that value for sigma
        if data_sig=="p":  self.sig = func.poisson_errors(self.y)
        elif type(data_sig)==type(""):  self.sig = eval(data_sig)*sc.ones(len(self.x))
        else:  self.sig = data_sig
        self.name = "fit"
        self.function = func.gaussian_function
        """ self.evaluator = R_squared """
        self.params = [1.0, 1.0, 1.0] # starting params - to be best params
        self.step = [0.1, 0.1, 0.1]  #fitting step sizes
        self.range = [None, None, None]  # (<#>, <min>, <max>, <step_size>), for surface
        self.min = [None, None, None]  # parameter min bounds
        self.max = [None, None, None]  # parameter max bounds
        self.errors = [None, None, None]  # Fit errors, to be determined
        self.RR = R_squared(self, self.params)  # best r-squared value
        # secondary values
        self.plot_type = "scatter"  #or "bar"
        self.plot_axes = ["x", "y", 0]
        self.param_names = ["parameter1", "parameter2", "parameter3", "parameter4"]
    def plot(self, new_params=None):
        """ Plots the data with current params or input params;  ERRORBAR? """
        if new_params==None:  new_params=self.params
        plt.figure(1)
        if self.plot_type=="scatter":
            plt.scatter(self.x,self.y,2,'k','o')
        elif self.plot_type=="bar":
            plt.bar(self.x, self.y,
                    ((np.ma.max(self.x)-np.ma.min(self.x))/len(self.x)), align='center')
        else:
            print "!!! '{0}' is not a valid plot type!".format(self.plot_type)
            sys.exit(2)
        plt.plot(self.x, self.function(self.x, new_params), 'k--')
        plt.xlabel(self.plot_axes[0])
        plt.ylabel(self.plot_axes[1])
        plt.title(str(new_params)+str(R_squared(self, new_params)), fontsize=8 )
        out_name = self.name+"_plot_"+str(self.plot_axes[2])+".ps"
        plt.show()
        #plt.savefig(out_name, papertype='letter')
        plt.close('all')
        # counter allows for a series of uniquely-named plots
        self.plot_axes[2] = self.plot_axes[2]+1
    def update_params(self, new_params):
        """ Ensures R-squared is updated with parameters """
        self.params = new_params
        self.RR = R_squared(self, self.params)

def R_squared(fitter, new_params):
    #return func.R_squared_plane(fitter.x, new_params)
    new_y = fitter.function(fitter.x, new_params)
    pieces = (fitter.y - new_y)*(fitter.y - new_y) / (fitter.sig*fitter.sig)
    return sc.sum(pieces)

def get_surface(fitter, p1=0, p2=1, plot=0, path=None):  
    """ Returns a grid of R-squared values by varying 2 parameters in the given ranges
        param1,2 are (<param#>, <min>, <max>, <step_size>)
        path=[p1, p2] will draw the search path over the surface"""
    param1 = fitter.range[p1]
    param2 = fitter.range[p2]
    steps1 = int((param1[2]-param1[1])/param1[3])
    steps2 = int((param2[2]-param2[1])/param2[3])
    surface = sc.zeros((steps1, steps2))
    # generate test parameter set
    new = []
    for i in range(len(fitter.params)):
        new.append(fitter.params[i])
    new[param1[0]] = param1[1]  #set to min
    new[param2[0]] = param2[1]  #set to min
    for i in range(steps1):
        for j in range(steps2):
            surface[i,j] = R_squared(fitter, new)  # 'flipped' y-axis
            new[param2[0]] = new[param2[0]] + param2[3]
        new[param1[0]] = new[param1[0]] + param1[3]
        new[param2[0]] = param2[1]  #reset to min
    if plot!=0:  plot_surface(fitter, surface, param1, param2, plot, path)
    return surface

def plot_surface(fitter, surface, param1, param2, plot=1, path=None):
        # get location of min
        min, x_min, y_min = 10000, 0, 0
        for i in range(surface.shape[0]):
            for j in range(surface.shape[1]):
                if surface[i,j] < min:  min, x_min, y_min = surface[i,j], j, i
        # the mathematical gynamstics are required to get things to line up on imshow
        plt.figure(1)
        plt.imshow(np.sqrt(surface), origin='lower')
        plt.colorbar()
        plt.scatter(x_min, y_min, 2,'k','o')
        if path != None:  plt.plot((path[1]-param2[1])/param2[3], (path[0]-param1[1])/param1[3], "k")        
        # set labels
        xlabels, ylabels = [], []
        xlocs = sc.arange(0.0, param2[2]-param2[1], (param2[2]-param2[1])/5.0 ) / param2[3]
        for i in range(len(xlocs)):  xlabels.append(str(xlocs[i]*param2[3]+param2[1])[:5])
        plt.xticks(xlocs, xlabels)
        ylocs = sc.arange(0.0, param1[2]-param1[1], (param1[2]-param1[1])/5.0 ) / param1[3]
        for i in range(len(ylocs)):  ylabels.append(str(ylocs[i]*param1[3]+param1[1])[:5])
        plt.yticks(ylocs, ylabels)
        #plt.figure(2)
        # for Contour plot overlay
        X = (sc.arange(param2[1], param2[2], param2[3]) / param2[3]) - (param2[1]/param2[3])
        Y = (sc.arange(param1[1], param1[2], param1[3]) / param1[3]) - (param1[1]/param1[3])
        if len(X) > surface.shape[1]:  X = X[:-1]  # Needed in case of rounding error
        if len(Y) > surface.shape[0]:  Y = Y[:-1]  # Needed in case of rounding error
        CS = plt.contour(X, Y, surface, levels=[0.0, 10.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0])
        plt.clabel(CS, inline=1, fontsize=10)
        plt.xlabel(fitter.param_names[param2[0]])
        plt.ylabel(fitter.param_names[param1[0]])
        if type(plot)==type(""):  plt.savefig(plot+".png", papertype="letter")
        else:   plt.show()
        plt.close('all')

def gradient_descent(fitter, line_search=0, f=0.0001, its=1000, precision=0.000001,
                     scale=1, verbose=0):
    """ Super version from Bevington & Robinson;
        f is the gradient sensitivity
        its is the max # of iterations"""
    validate = 0  # check to see if validation is needed
    for i in range(len(fitter.min)):
        if fitter.min[i] != None:  validate = 1
    for i in range(len(fitter.max)):
        if fitter.max[i] != None:  validate = 1
    path, new = [], [] # Initialize new, path
    for i in range(len(fitter.params)):  new.append(fitter.params[i])
    for i in range(len(fitter.params)):  app_path(path, fitter.params)
    l = len(new)
    # get gradient
    # scale length, counter, move counter, difference in RRs
    S, loop, moves, diff, overrun = 10.0, 0, 0, 1000.0, 0  
    while ((abs(diff) > precision) or (loop == 0)):
        denom, partial = 0.0, []
        for i in range(l):
            new[i] = new[i] + f*fitter.step[i]     # iterate parameter
            RR_grad = R_squared(fitter, new)    # get new R-squared
            partial.append((RR_grad-fitter.RR)/f)  # append the partial derivative
            denom = denom + partial[i]*partial[i]  # build the normalizing factor
            new[i] = fitter.params[i]           # reset the changed parameter
        denom = np.sqrt(denom)
        for i in range(l):  new[i] = new[i] -1.0*S*(partial[i]/denom)*fitter.step[i]  # EQN 8.20, p.154
        if validate == 1:   # Check new params against bounds
            if is_move_valid(fitter, new, verbose) == 0:
                S = S*0.5   # The only way to change S when scale==0
                continue
        RR_new = R_squared(fitter, new)
        diff = RR_new - fitter.RR
        if diff < 0.0:      # Better Fit
            if line_search==1:  new, RR_new = line_search(fitter, new, RR_new, validate, verbose)
            if scale==1:  S = S*1.2     # Increase scale factor to move faster
            if verbose==1:  print "# - Moving from \n{0} to \n{1}, R-squared={2}, diff of {3}".format(fitter.params, new, RR_new, (RR_new - fitter.RR))
            app_path(path, new)
            for i in range(l):  fitter.params[i] = new[i]
            fitter.RR = RR_new
            moves = moves + 1
        else:               # Worse fit
            if scale==1:  S = S*0.8     # Decrease scale factor to move less
            if verbose==1:  print "\n##- No Move from \n{0} to \n{1}, diff of {2}, {3}".format(fitter.params, new, diff, S)
            for i in range(l):  new[i] = fitter.params[i]  # reset new
        # Iterate and check for overrun
        loop = loop + 1
        if (loop > its):
            print "!!! - Exited by passing iteration threshold: {0}".format(loop)
            overrun = 1
            break
    if overrun == 0:
        print "# - Gradient Descent successful after {0} iterations, {1} moves".format(loop, moves)
    print "# - Final values:  {0} : {1}".format(fitter.params, fitter.RR)
    return sc.array(path)

def MCMC(fitter, iters=1000, annealing=0, verbose=0, seed=None):
    """ Monte Carlo - Markov Chain, with Metropolis-Hastings step size
        iters is maximum number of itereations, dropout is success ratio threshold
        to end search;  anealing turns on simulated anealing method
        Step sizes of 1/20 - 1/40 of the order of the expected parameters are good."""
    if seed == None:  np.random.seed(seed)
    else:  np.random.seed(time.time())
    validate = 0  # check to see if validation is needed
    for i in range(len(fitter.min)):
        if fitter.min[i] != None:  validate = 1
    for i in range(len(fitter.max)):
        if fitter.max[i] != None:  validate = 1
    # Initialize new params, path holder, loop counter
    path, new, p_best = [], [], []
    loops, inval, moves = 0, 0, 0
    for i in range(len(fitter.params)):  new.append(fitter.params[i])
    for i in range(len(fitter.params)):  p_best.append(fitter.params[i])
    for i in range(len(fitter.params)):  app_path(path, fitter.params)
    RR_best = R_squared(fitter, fitter.params)
    l = len(new)
    while loops < iters:
        for i in range(l):
            new[i] = np.random.normal(fitter.params[i], fitter.step[i])
        if validate == 1:   # Check new params against bounds, get new params if bad
            if is_move_valid(fitter, new, verbose) == 0:
                inval = inval + 1
                if inval > 100:  print "!!! - unable to find valid step.  Exiting..."; break
                continue
        inval = 0
        RR_new = R_squared(fitter, new)
        #print RR_new
        if RR_new < RR_best:
            for i in range(l):  p_best[i] = new[i]
            RR_best = R_squared(fitter, new)
        quot = np.exp(-(RR_new - fitter.RR)/ 2.0)  # quot throws nans!
        #print quot, RR_new, fitter.RR
        #if ma.isnan(quot):  print quot, fitter.RR, np.exp(-fitter.RR/2.0), RR_new, np.exp(-RR_new/2.0)
        rand = np.random.uniform()
        #print rand, 2.0*np.log(rand)
        if quot > rand:  # accept point
            if verbose==1:  print "# - Moving from \n{0} to \n{1}, R-squared={2}, diff of {3}".format(fitter.params, new, RR_new, (RR_new - fitter.RR))
            for i in range(l):  fitter.params[i] = new[i]
            fitter.RR = RR_new
            moves = moves + 1
        else:               # Worse fit
            if verbose==1:  print "\n##- No Move from \n{0} to \n{1}, diff of {2}".format(fitter.params, new, (RR_new - fitter.RR))
        app_path(path, fitter.params)
        # Iterate and check for overrun
        loops = loops + 1
    #if (loops > 99) and (float(moves)/float(loop) < dropout):  #IMPLEMENT DROPOUT?
    #    print "# Exiting due to success dropout condition"
    # Print best params
    print "# - MCMC exit with best params, R-squared: {0}; {1}".format(p_best, RR_best)
    # get params from posterior distribution
    path = sc.array(path)
    p_mean, p_std = [], []
    for i in range(l):  p_mean.append(sc.mean(path[:,i])); p_std.append(sc.std(path[:,i]))
    print "# - After {0} moves, {1} iterations".format(moves, iters)
    print "# - Final values (mean):  {0}".format(p_mean)
    print "# - Value deviations (stdev):  {0}".format(p_std)
    # which parameter set should be returned as best?
    return path

def app_path(path, params):
    holder=[]
    for i in range(len(params)):  holder.append(params[i])
    path.append(holder)

def is_move_valid(fitter, new, verbose=0):
    for i in range(len(fitter.params)):
        if fitter.min[i] != None:
            if new[i] < fitter.min[i]:
                if verbose==1:
                    print "! - parameter {0} is below bound! ({1})".format(i, new[i])
                return 0
        if fitter.max[i] != None:
            if new[i] > fitter.max[i]:
                if verbose==1:
                    print "! - parameter {0} is above bound! ({1})".format(i, new[i])
                return 0
    return 1

def line_search(fitter, new, RR_new, validate=0, verbose=0):
    n_params = len(new)
    params = [fitter.params, new]
    diff = RR_new - fitter.RR
    tries = 0
    while diff < 0.0:
        new_params=[]  # new params are twice as far as last set:  X" = 2x'- x
        for i in range(n_params):  new_params.append(2.0*params[-1][i] - params[0][i])
        if validate==1:
            if is_move_valid(fitter, new_params, verbose=0)==0:  break
        RR_newer = R_squared(fitter, new_params)
        diff = RR_newer - RR_new
        print new_params, diff
        if diff < 0.0:  # better fit than last
            params.append(new_params)
            RR_new = RR_newer
        tries = tries + 1
        if tries > 100:  print "! - Too many steps in line search";  break
    print "\n"
    return params[-1], RR_new


# conjugate gradient descent (http://en.wikipedia.org/wiki/Conjugate_gradient_method)

def get_Hessian_errors(fitter, like=0):
    """ From Hessian script """
    errors = he.get_hessian_errors(fitter.function, fitter.params, fitter.x, fitter.y,
                          fitter.sig, fitter.step, like)
    return errors

def get_correlation(x, y):
    """ returns correlation 'r' from x and y values.
        see Bevington, 'Data Reduction and Error Analysis', EQN 11.17"""
    N = len(x)
    sum_x, sum_y = sc.sum(x), sc.sum(y)
    sum_xy, sum_xx, sum_yy = sc.sum(x*y), sc.sum(x*x), sc.sum(y*y)
    top = N*sum_xy - (sum_x*sum_y)
    botx = sc.sqrt(N*sum_xx - (sum_x*sum_x) )
    boty = sc.sqrt(N*sum_yy - (sum_y*sum_y) )
    return top/(botx*boty)

# reduced chi-squared

if __name__ == "__main__":
    data = fi.read_data("Pal13_backsub_cluscut_cut.txt")#"eps_garb.txt")
    fitter=ToFit(data, cols=[1,0,"p"])
    fitter.function=func.get_2gauss_y
    fitter.update_params([4.0, 0.3, 1.0, 8.0])
    fitter.step = [0.01, 0.01, 0.01, 0.01]
    #fitter.range = [(0, -2.0, 2.0, 0.1), (1, 0.0, 0.2, 0.005), (2, -10.0, 0.0, 0.25)]
    fitter.range = [(0, 4.0, 5.0, 0.05), (1, 0.2, 1.2, 0.05), (2, 0.2, 1.2, 0.05), (3, 10.0, 20.0, 0.5)]
    print fitter.params, fitter.RR
    path = gradient_descent(fitter, its=10000, line_search=0)
    fitter.plot()
    surface = get_surface(fitter, 0, 1, plot=1, path=(path[:,0], path[:,1]))
    #surface = get_surface(fitter, (1, 1.0, 6.0, 0.1), (2, 0.5, 1.1, 0.01),
    #                      plot=1, path=(path[:,1], path[:,2]))
    print "# --- Done"
