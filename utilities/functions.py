import math as m
import numpy as np
import scipy as sc

'''python script for containing useful functions, esp. for curve fitting.
See "make_data.py" for creation of test data sets.
Matthew Newby, May 2, 2011'''

# Remake functions as class types?  defs: "generate_y, derivative...?"        

def gaussian_function(x_in, params):
    # 3 params, fits the amplitude SQRT!!!, mu, and stdev to a symmetric gaussian
    amp, mu, sigma = params
    y_out = sc.zeros(len(x_in))
    sig_sq = sigma*sigma
    #factor = 1.0 /(np.sqrt(2.0*m.pi*sig_sq))
    for i in range(len(x_in)):
        exponent = -0.5*( ((x_in[i] - mu)*(x_in[i] - mu)) / sig_sq )
        y_out[i] = amp*amp*np.exp(exponent)
    return y_out

def gauss_plus_floor(x_in, params):
    # 3 params, fits the amplitude SQRT!!!, mu, and stdev to a symmetric gaussian
    amp, mu, sigma, floor = params
    y_out = sc.zeros(len(x_in))
    sig_sq = sigma*sigma
    #factor = 1.0 /(np.sqrt(2.0*m.pi*sig_sq))
    for i in range(len(x_in)):
        exponent = -0.5*( ((x_in[i] - mu)*(x_in[i] - mu)) / sig_sq )
        y_out[i] = amp*amp*np.exp(exponent) + floor
    return y_out

def get_2gauss_y (x, parameters):
    #Creates Y values for given x values for a double-sided gaussian function.
    # that is, a Gaussian with different sigma's on either side of the mean
    # From "gaus_fit_enc.py"
    length = len(x)
    mu, sigma_left, sigma_right, amplitude = parameters
    line_y = sc.zeros(length)
    #Create the gaussian from parameter sets
    for i in range(length):
        if x[i] < mu:  sigma = sigma_left
        else:  sigma = sigma_right
        exponent = -0.5 * ( (x[i] - mu) / sigma )**2
        stuff = m.exp(exponent)
        line_y[i] = stuff*amplitude
    return line_y

def two_stage (x, params):
    y = sc.zeros(len(x))
    for i in range(len(x)):
        first = (1.0/(1.0 + np.exp( (x[i]-params[1])/(params[2]) )))
        if x[i] > params[3]:
            second = np.sqrt((x[i]-params[3]) / params[4])
        else:  second = 0.0
        #    second = y[i] * ( 1.0 / (1.0 + ( (x[i]-params[3]) / params[4]) ) )
        #else:  second = 1.0
        y[i] = params[0]*first - second
    return y

def plane_fit_old(data, params):
    """ Attempting to output distance 'd' as y, x=array([x,y,z]) fit points.
        y = an array of 0's.  Might have to circumvent sigma as well...
        Looks like the plotting function is the only thing that will be affected"""
    # 4 parameters:  a, b, c, d gives plane:  ax +by +cz +d = 0
    a, b, c, d = params
    x, y, z = data[:,0], data[:,1], data[:,2]
    if len(data[0,:]) > 3:  errors= data[:,3]
    else:  errors = sc.ones(len(data[:,0]) )
    # Plane normal:  V = [a,b,c]
    bottom = np.sqrt(a*a + b*b + c*c)
    D = sc.zeros(len(x))
    # Distance, Wolfram Mathworld
    for i in range(len(x)):
        top = a*x[i] + b*y[i] + c*z[i] + d
        D[i] = top / (bottom*errors[i])
    return D

def plane_fit(data, params):
    """ As above, except that parameters are theta, phi, d
        angles are postion of normalized plane normal on unit
        sphere, in radians"""
    theta, phi, d = params
    a, b, c = sc.cos(theta)*sc.sin(phi), sc.sin(theta)*sc.sin(phi), sc.cos(phi)
    x, y, z = data[:,0], data[:,1], data[:,2]
    D = sc.zeros(len(x))
    for i in range(len(x)):  D[i] = a*x[i] + b*y[i] + c*z[i] + d
    return D
    
def plummer_density(x, params):
    stuff = (1.0 + (x*x)/(params[1]*params[1]) )**(-5/2)
    return (params[0])/(params[1]*params[1]*params[1])*stuff

    
def double_gaussian(x, params):
    # 6 parameters
    A1, mu1, sig1, A2, mu2, sig2 = params
    G1 = A1*A1*sc.exp(-1.0*(x-mu1)*(x-mu1) / (2.0*sig1*sig1) )
    G2 = A2*A2*sc.exp(-1.0*(x-mu2)*(x-mu2) / (2.0*sig2*sig2) )
    #y = (0.315*x) + 34.2 #offset for virgo, 20-20.5 cut
    #y = (2.0*x) + 160.0 #offset for virgo, full cut
    return G1 + G2 #+ y


def double_gauss_line(x, params):
    # 8 parameters
    A1, mu1, sig1, A2, mu2, sig2, aa, bb = params
    G1 = A1*A1*sc.exp(-1.0*(x-mu1)*(x-mu1) / (2.0*sig1*sig1) )
    G2 = A2*A2*sc.exp(-1.0*(x-mu2)*(x-mu2) / (2.0*sig2*sig2) )
    y = (aa*x) + bb
    return G1 + G2 + y


def double_gaussian_one_fixed(x, params):
    # 6 parameters, 4 fit;  holds the mean and sigma constant for 2nd gaussian
    A1, mu1, sig1, A2, mu2, sig2 = params
    G1 = A1*A1*sc.exp(-1.0*(x-mu1)*(x-mu1) / (2.0*sig1*sig1) )
    G2 = A2*A2*sc.exp(-1.0*(x-mu2)*(x-mu2) / (2.0*sig2*sig2) )
    return G1 + G2 


def quad_fat_gaussians(x, params):
    """ Fits two pairs of Gaussians to data;  each pair has the same mean.
        10 parameters """
    A1, mu1, sig1, A2, mu2, sig2, A3, sig3, A4, sig4 = params
    G1 = A1*A1*sc.exp(-1.0*(x-mu1)*(x-mu1) / (2.0*sig1*sig1) )
    G2 = A2*A2*sc.exp(-1.0*(x-mu2)*(x-mu2) / (2.0*sig2*sig2) )
    G3 = A3*A3*sc.exp(-1.0*(x-mu1)*(x-mu1) / (2.0*sig3*sig3) )
    G4 = A4*A4*sc.exp(-1.0*(x-mu2)*(x-mu2) / (2.0*sig4*sig4) )
    #y = (0.315*x) + 34.2 #offset for virgo, 20-20.5 cut
    #y = (2.0*x) + 160.0 #offset for virgo, full cut
    return G1 + G2 + G3 + G4 #+ y


def quad_fat_gauss_line(x, params):
    """ Fits two pairs of Gaussians to data;  each pair has the same mean.
        12 parameters """
    A1, mu1, sig1, A2, mu2, sig2, A3, sig3, A4, sig4, aa, bb = params
    G1 = A1*A1*sc.exp(-1.0*(x-mu1)*(x-mu1) / (2.0*sig1*sig1) )
    G2 = A2*A2*sc.exp(-1.0*(x-mu2)*(x-mu2) / (2.0*sig2*sig2) )
    G3 = A3*A3*sc.exp(-1.0*(x-mu1)*(x-mu1) / (2.0*sig3*sig3) )
    G4 = A4*A4*sc.exp(-1.0*(x-mu2)*(x-mu2) / (2.0*sig4*sig4) )
    y = (aa*x) + bb
    return G1 + G2 + G3 + G4 + y


def gaussian_constmean(x_in, params):
    """Fits only sigma and amplitude, holding the mean constant"""
    mu = 0.0
    amp, sigma = params
    y_out = sc.zeros(len(x_in))
    sig_sq = sigma*sigma
    #factor = 1.0 /(np.sqrt(2.0*m.pi*sig_sq))
    for i in range(len(x_in)):
        exponent = -0.5*( ((x_in[i] - mu)*(x_in[i] - mu)) / sig_sq )
        y_out[i] = amp*amp*np.exp(exponent)
    return y_out


def double_sigmoid(x, params):
    y = sc.zeros(len(x))
    for i in range(len(x)):
        first = (1.0/(1.0 + np.exp( (x[i]-params[1])/(params[2]) )))
        second = (1.0/(1.0 + np.exp( (x[i]-params[3])/(params[4]) )))
        y[i] = params[0]*first*second
    return y

def inverse_stable(x, params):
    y = 1.0 / (1.0 + ( (x-params[1]) / params[0]) ) 
    return y

def Nth_order_polynomial(x, params):
    #fits a polynomial of order=number of parameters - 1
    y = sc.zeros(len(x))
    for i in range(len(x)):
        y[i] = params[0]
        for N in range(1,len(params)):
            holder = 1.0
            for j in range(N):
                holder = holder*x[i]
            y[i] = y[i] + params[N]*holder
    return y

def plank_curve(x, params):
    y = sc.zeros(len(x))
    for i in range(len(x)):
        x_new = x[i] - params[2]
        holder = 1.0 / (np.exp(params[1]/x_new) - 1.0)
        y[i] = params[0]*holder/(x_new*x_new*x_new)
    return y

def critical_dampning(x,params):
    # 3 parameters
    y = sc.zeros(len(x), float)
    for i in range(len(x)):
        y[i] = (params[0] + x[i]*params[1])*np.exp(-1.0*params[2]*x[i])
    return y

def synaptic_response(x, params):
    #http://rsif.royalsocietypublishing.org/content/5/29/1429.full
    y = sc.zeros(len(x), float)
    const = np.exp(-1.0*params[0])
    for i in range(len(x)):
        front = params[0]*params[0]*np.exp(-1.0*params[0]*x[i]) / (1.0-const)
        back = x[i] + (const / (1.0-const))
        y[i] = front * back
    return y

def average_fit(x, params):
    #One parameter, fits an average to error-weighted data
    y = params[0]*sc.ones(len(x))
    return y

def sigmoid_function(x, params):
    #Looks good, takes three parameters
    data = ( params[0]/(1 + np.exp(-1.0*(x-params[1])) ) ) + params[2]
    return data

def reversed_sigmoid_function(x, params):
    #Looks good, takes three parameters; exp argument not multiplied by -1
    data = ( params[0]/(1 + np.exp(1.0*(x-params[1])) ) ) + params[2]
    return data

def SDSS_detection_efficiency(x):
    s = [0.9402, 1.6171, 23.5877]
    detection_efficiency = s[0] / (np.exp(s[1]*(x - s[2])) + 1.)
    return detection_efficiency

def mod_sigmoid_function(x, params):
    #Looks good, takes three parameters
    data = ( params[0]/(1 + np.exp(-1.0*(x-params[1])) ) ) + params[2] + (x*params[3])
    return data

def parabola(x, params):
    #takes 2 params
    data = 1.0 + params[0]*((x - params[1])**2)
    return data

def linear_function(x, params):
    #looks like this function is ok
    y = (x*params[0] + params[1])
    return y

def const_linear_function(x, params):
    #looks like this function is ok
    #y = ((x-20.0)*params[0] + 1.24127)
    y = ((x-18.0)*params[0] + 1.24009)
    return y

def make_line_errors(m, b, yerr=0.0, xerr=0.0, xrange=[0.0,100.0], step=1.0):
        np.random.seed( int(time.time()) )
        length = int((xrange[1]-xrange[0])/step)
        data = sc.zeros((length,4))
        for x in range(length):
                y = m*x + b
                if (xerr <= 0):  data[x,0] = x
                else:  data[x,0] = np.random.normal(x, xerr)
                if (yerr <= 0):  data[x,1] = y
                else:  data[x,1] = np.random.normal(y, yerr)
        data[:,2] = xerr
        data[:,3] = yerr
        return data

def gaussian_function_normal(x_in, params):
    #this function appears to be crap when used to fit, but produces the right profile otherwise...
    y_out = sc.zeros(len(x_in))
    sig_sq = params[1]*params[1]
    factor = 1.0 /(np.sqrt(2.0*m.pi*sig_sq))
    for i in range(len(x_in)):
        exponent = -0.5*( ((x_in[i] - params[0])*(x_in[i] - params[0])) / sig_sq )
        y_out[i] = factor*np.exp(exponent)
    return y_out

def make_gauss_errors(mu, sigma, yerr=0.0, xerr=0.0, xrange=[0.0, 10.0], step=0.1):
    np.random.seed( int(time.time()) )
    length = int((xrange[1]-xrange[0])/step)
    data = sc.zeros((length,4))
    x = sc.arange(xrange[0], xrange[1], step)
    y = gaussian_function(x, [mu, sigma])
    for i in range(length):
        if (xerr <= 0):  data[i,0] = x[i]
        else:  data[i,0] = np.random.normal(x[i], xerr)
        if (yerr <= 0):  data[i,1] = y[i]
        else:  data[i,1] = np.random.normal(y[i], yerr)
    data[:,2] = xerr
    data[:,3] = yerr
    return data

def SDSS_error_function(x_in, params):
    if len(params) != 3:  print 'wrong number of input parameters!'
    exponent = (x_in*params[1] + params[2])
    y = params[0] + np.exp(exponent)
    return y

def SDSS_error_function2(x_in, params):
    if len(params) != 2:  print 'wrong number of input parameters!'
    y = (x_in*params[0] + params[1])
    return y

def poisson_errors(n, rigorous=False):
    #upper error bar
    #lambda_up = n + np.sqrt(n+1.0) + 1.0
    #lower error bar
    #lambda_down = n*( (1.0 - (1.0/(9.0*n)) - (1.0/(3.0*np.sqrt(n))))**3 )
    #poisson error:
    lambda_up = np.sqrt(n+1.0) + 1.0
    return lambda_up

def integrate_gaussian(x1, x2, mu=0.0, sig=1.0, amp=1.0):
    """ Provides the area under a gaussian between two points. 
        Defaults to a standard normal distribution"""
    import scipy.special as ss
    sq2 = m.sqrt(2)
    sqpi = m.sqrt(m.pi)
    n1 = (x1 - mu)/sig
    n2 = (x2 - mu)/sig
    prefix = 0.5*amp*abs(sig)*sqpi 
    return prefix*(ss.erf(n2/sq2) - ss.erf(n1/sq2))
