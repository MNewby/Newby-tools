#! /usr/bin/env python  #'Bang' line - modify as needed

import sys
sys.path.insert(0, '../utilities')
import time
import math as ma
import numpy as np
import scipy as sc
import files as fi
import progress as pr
import astro_coordinates as co
#from __future__ import print_function


"""This script creates test stripes for stream mle code and milkyway@home.
Matthew Newby (RPI), March 28, 2011
"""
deg = 180.0 / ma.pi 
rad = ma.pi / 180.0 

# Finish stream_length, convolve, and master function

class ParamSet:
    """ A set of milkyway@home fit parameters"""
    def __init__(self, paramString=None):
    # SDSS Stripe Number
        self.wedge = 82
    # Background weight (epsilon), 4 parameters (?, q, r_0, ?)
        self.background = [0.0, 1.0, 0.455, 19.5, 1.0]
    # stream weight (epsilon), 5 parameters (mu, r, theta, phi, sigma) per stream
        self.streams = [] #[ [-2.0, 360.0, 21.0, 0.0, 0.0, 1.0] ]#,
                #[-4.0, 330.0, 10.0, 18.0, 0.0, 1.0],
                #[-10.0, 400.0, 10.0, 20.0, 0.0, 5.0]
                #]
    # stripe parameters:  mu, nu, r; min, max, n_steps <-n_steps doesn't seem to be used...
        self.stripe = [(310.0, 419.0, 10),(-1.25, 1.25, 10),(16.0,22.5,10)]
        if paramString != None:
            self.load_params(paramString)
        self.update_refs()
    def update_refs(self):
    # Easy external references for parameters
        self.back_weight = self.background[0]
        self.q = self.background[2]
        self.r0 = self.background[3]
        self.num_streams=len(self.streams)
        if self.num_streams > 0:
            self.stream_weight = [self.streams[0][0]]
            self.mu = [self.streams[0][1]]
            self.R = [self.streams[0][2]]
            self.theta = [self.streams[0][3]]  # IN DEGREES!!!!
            self.phi = [self.streams[0][4]]   # IN DEGREES!!!!
            self.sigma = [self.streams[0][5]]
            if len(self.streams) > 1:
                for i in range(1,len(self.streams)):
                    self.stream_weight.append(self.streams[i][0])
                    self.mu.append(self.streams[i][1])
                    self.R.append(self.streams[i][2])
                    self.theta.append(self.streams[i][3])
                    self.phi.append(self.streams[i][4])
                    self.sigma.append(self.streams[i][5])
        self.mu_lim = self.stripe[0]
        self.nu_lim = self.stripe[1]
        # Damn r in param files is actually g, will convert here to make referencing easy:
        self.r_lim = (co.getr(self.stripe[2][0]), co.getr(self.stripe[2][1]), self.stripe[2][2])
    def load_params(self, params):
        self.background = params[0:4]
        self.streams = []
        for i in range(params[5:]/6 ):
            holder = params[(i+5):(i+11)]
            self.streams.append(holder)
    def print_params(self):
        print "\n!!! NEEDS UPDATING TO CONFORM WITH PARAM FILE FORMATS !!!\n"
        print "background_weight:", self.background[0]
        print "background_parameters[4]:", self.background[1:]
        for i in range(len(self.streams)):
            print "stream_weight:", self.streams[i][0]
            print "stream_parameters[5]:", self.streams[i][1:]
        print self.stripe


def hernquist_profile(x,y,z, q,r0):
    """ Galactic-centric hernquist density profile """
    r = sc.sqrt( (x*x) + (y*y) + ((z*z)/(q*q)) )
    inv_rho = r*(r + r0)*(r + r0)*(r + r0)
    return (1.0/inv_rho)

def generate_stream(num_stars, u_min, u_max, sigma):
    """ stream generator, uses active generation technique """
    #low, high = -1.0*(length/2.0), (length/2.0)
    u = np.random.normal(0.0, sigma, size=num_stars)
    v = np.random.normal(0.0, sigma, size=num_stars)
    w = np.random.uniform(u_min, u_max, num_stars)
    return u,v,w

def stream_into_stripe(params, sn, N_stars, batch=1000, fileout="streamgen82.txt",
                       detection=1, convolve=1, append=0, progbar=1):
    """ sn is stream number"""
    # Initialize file
    if append==0:
        out = open(fileout, 'w')
        out.write("# "+str(N_stars)+" stars, l,b,r \n")
        out.close()
    # Get constants
    mu,R,theta,phi,sigma,wedge = \
        params.mu[sn],params.R[sn],params.theta[sn],params.phi[sn],params.sigma[sn], params.wedge
    u_min, u_max = get_stream_length(params, sn, accuracy=0.0001)
    nu_min, nu_max = params.nu_lim[0], params.nu_lim[1]
    mu_min, mu_max = params.mu_lim[0], params.mu_lim[1]
    g_min, g_max = params.stripe[2][0], params.stripe[2][1]
    #print "# - Generating Stream {0}, using parameters {1}, {2}, {3}, {4}, {5}".format(
    #    sn, mu, R, theta, phi, sigma)
    N_out = 0
    pb = pr.Progressbar(steps=N_stars, prefix="Stream {0} progress:".format(sn), 
        suffix="Generating {0} Stars".format(N_stars), symbol="#", 
                active="=", brackets="[]", percent=True, size=40)
    while N_out < N_stars:
        mu_killed, nu_killed, mu_saved = 0,0,0
        u,v,w = generate_stream(batch, u_min, u_max, sigma)
        holder = []
        for i in range(len(u)):
            mu1, nu1, r1 = co.streamToGC(u[i],v[i],w[i],mu,R,theta*deg,phi*deg,wedge)
            if (nu1 < nu_min) or (nu1 > nu_max):  nu_killed=nu_killed+1; continue
            if (mu_max > 360.0):
                if (mu1 > (mu_max-360.0)) and (mu1 < mu_min):  mu_killed=mu_killed+1; continue
            else:
                if (mu1 < mu_min) or (mu1 > mu_max):  mu_killed=mu_killed+1; continue
            # Convolve
            if convolve==1:
                r1 = star_convolution(r1)
            # Detection
            if detection==1:
                m_g = co.getg(r1)
                if np.random.uniform() > (sigmoid_error(m_g)):  continue
            if (co.getg(r1) < g_min) or (co.getg(r1) > g_max): continue
            # When a point passes all the testsm add it to the set
            l,b,r1 = co.GC2lbr(mu1, nu1, r1, wedge)
            #co.stream2xyz(u[i],v[i],w[i],mu,R,theta,phi,wedge)
            holder.append([round(l,6),round(b,6),round(r1,6)])
            N_out = N_out + 1
        if N_out > N_stars:
            slice = -1*(N_out-N_stars)
            holder = holder[:slice]
            #print "#---Sliced {0} stars to make quota".format(str(-1*slice))
            N_out = N_out + slice #Slice is negative
        #append to file
        if len(holder) != 0:
            if fi.append_data(sc.array(holder), fileout, delimiter=" ") == 1:
                #print "#---Stream Progress: {0} stars of {1} total stars generated".format(N_out, N_stars)
                #print "# !!! - mu killed: {0}, mu_saved: {1}, nu_killed: {2}".format(mu_killed, mu_saved, nu_killed)
                pb.updatebar(float(N_out)/float(N_stars))
    print "#---Stream {0} generation succeeded, written as {1}".format(sn, fileout)
    return fileout

def get_stream_length(params, N=0, accuracy=0.0001):
    """ gets the length of a stream within an SDSS wedge; N=stream number"""
    mu, R, theta, phi, wedge = params.mu[N],params.R[N],params.theta[N],params.phi[N],params.wedge
    # make a point (or points) in xyz along stream direction
    u1,v1,w1 = 0.1,  0.0, 0.0
    u2,v2,w2 = -0.1,  0.0, 0.0
    #check to see which is closest to the nu = +2.5 boundry, flip if necessary
    mu1, nu1, r1 = co.streamToGC(u1,v1,w1,mu, R, theta*deg, phi*deg, wedge)
    mu2, nu2, r2 = co.streamToGC(u2,v2,w2,mu, R, theta*deg, phi*deg, wedge)
    if np.fabs(nu1 - 2.5) > np.fabs(nu2 - 2.5):
        u1, u2 = (-1.0*u1), (-1.0*u2)
        temp = nu1
        nu1, nu2 = nu2, temp
    # check them against wedge boundries (mu, nu, r)
    test = 0
    while test==0:
        mu1, nu1, r1 = co.streamToGC(u1,v1,w1,mu, R, theta*deg, phi*deg, wedge)
        # account for wrap-around
        if params.mu_lim[1] > 360.0:
            if mu1 < params.mu_lim[0]:  mu1 = mu1 + 360.0
        if (mu1 < params.mu_lim[0]) or (mu1 > params.mu_lim[1]):
            #print "u1 ({0}) finished due to mu lim: {1}, {2}, {3}".format(u1, mu1, nu1, r1)
            break
        if (nu1 > 2.5):
            #print "u1 ({0}) finished due to nu lim: {1}, {2}, {3}".format(u1, mu1, nu1, r1)
            break
        if (nu1 < -2.5):
            #print "!!! WARNING:  u1 EXITED DUE TO OPPOSITE THRESHOLD"
            #print "u1 ({0}) finished due to nu lim: {1}, {2}, {3}".format(u1, mu1, nu1, r1)
            break
        if (r1 < params.r_lim[0]) or (r1 > params.r_lim[1]):
            #print "u1 ({0}) finished due to r lim: {1}, {2}, {3}".format(u1, mu1, nu1, r1)
            break
        u1 = u1 + 0.1
        if u1 > 100.0:  print "!!! u1 {0} reached threshold limit!".format(u1); test=1
    # Do for u2
    test = 0
    while test==0:
        mu2, nu2, r2 = co.streamToGC(u2,v2,w2,mu, R, theta*deg, phi*deg, wedge)
        # account for wrap-around
        if params.mu_lim[1] > 360.0:
            if mu2 < params.mu_lim[0]:  mu2 = mu2 + 360.0
        if (mu2 < params.mu_lim[0]) or (mu2 > params.mu_lim[1]):
            #print "u2 ({0}) finished due to mu lim: {1}, {2}, {3}".format(u2, mu2, nu2, r2)
            break
        if (nu2 < -2.5):
            #print "u2 ({0}) finished due to nu lim: {1}, {2}, {3}".format(u2, mu2, nu2, r2)
            break
        if (nu2 > 2.5):
            #print "!!! WARNING:  u2 EXITED DUE TO OPPOSITE THRESHOLD"
            #print "u2 ({0}) finished due to nu lim: {1}, {2}, {3}".format(u2, mu2, nu2, r2)
            break
        if (r2 < params.r_lim[0]) or (r2 > params.r_lim[1]):
            #print "u2 ({0}) finished due to r lim: {1}, {2}, {3}".format(u2, mu2, nu2, r2)
            break
        u2 = u2 - 0.1
        if u2 < -100.0:  print "!!! u1 {0} reached threshold limit!".format(u1); test=1
    # finish up
    length = np.fabs(u1-u2)
    print "# Stream is {0} kpc long within wedge boundaries".format(length)
    while np.fabs(u1-u2) < 20.0:
        u1 = u1*2.0
        u2 = u2*2.0
    if u2 > u1:  return u1, u2
    else:        return u2, u1

# Need to cut at borders?, apply convolution, rolloff

def generate_background(num_stars, params, batch=1000, fail_quit=100,
                        fileout="Backgen82.txt", detection=1, convolve=1, append=0):
    """Density of smooth halo background, as a function of position"""
    # Initialize file
    if append==0:
        out = open(fileout, 'w')
        out.write("# "+str(num_stars)+" stars, l,b,r \n")
        out.close()
    # Get integral
    tot_prob = get_max_prob(params)  #rough3d_integrator(params)
    N_out, fails = 0, 0
    g_min, g_max = params.stripe[2][0], params.stripe[2][1]
    while N_out < num_stars:
        # Generate Points
        mu, nu, r = get_stripe_points(params.mu_lim, params.nu_lim, params.r_lim, batch)
        # Test points for inclusion
        x,y,z = co.GC2xyz(mu, nu, r, params.wedge)
        holder = []
        for i in range(len(mu)):
            rho = hernquist_profile(x[i],y[i],z[i], params.q, params.r0)
            #print (rho / tot_prob), rho, tot_prob
            if (rho / tot_prob) > np.random.uniform():  
                l,b,r1 = co.xyz2lbr (x[i], y[i], z[i])
                # Convolve
                if convolve==1:
                    r1 = star_convolution(r1)
                # Detection
                if detection==1:
                    m_g = co.getg(r1)
                    if np.random.uniform() > (sigmoid_error(m_g)):  continue
                if (co.getg(r1) < g_min) or (co.getg(r1) > g_max): continue
                # Add to keepers
                holder.append([round(l,6),round(b,6),round(r1,6)])
                N_out = N_out + 1
        # Failure code
        if len(holder) == 0:  fails = fails + 1
        if fails >= fail_quit:  break
        # Remove possible excess stars
        if N_out > num_stars:
            slice = -1*(N_out-num_stars)
            holder = holder[:slice]
            print "#---Sliced {0} stars to make quota".format(str(-1*slice))
            N_out = N_out + slice #Slice is negative
        # Add to dataset
        if len(holder) != 0:
            if fi.append_data(sc.array(holder), fileout, delimiter=" ") == 1:
                print "#---Background Progress: {0} stars of {1} total stars generated".format(N_out, num_stars)
    if fails >= fail_quit:  
        print "!!! Background generation FAILED due to overstepping empty batch limit: {0}".format(fail_quit)
    else:
        print "#---Background generation succeeded, written as {0}, with {1} empty batches".format(fileout, fails)
    return fileout

def get_stripe_points(mu_lim, nu_lim, r_lim, number=1):   #GC checked
    mu = np.random.uniform(low=mu_lim[0], high=mu_lim[1], size=number)
    u, w = np.random.uniform(size=number), np.random.uniform(size=number)
    nu = sc.arcsin( u*(sc.sin(nu_lim[1])-sc.sin(nu_lim[0])) + sc.sin(nu_lim[0]) )
    r = r_lim[1]*(w**(1./3.))
    return (mu,nu,r)

def get_max_prob(params):
    mu_min, mu_max, mu_steps = params.mu_lim
    nu_min, nu_max, nu_steps = params.nu_lim
    r_min, r_max, r_steps = params.r_lim
    total_steps = ((mu_steps-1)*(nu_steps-1)*(r_steps-1))
    Dmu = (mu_max-mu_min) / mu_steps
    Dnu = (nu_max-nu_min) / nu_steps
    Dr = (r_max-r_min) / r_steps
    print "Integrating function,", total_steps, "steps"
    max_rho, count, progress = 0.0, 0, 0
    # Integrate from value at center of each bin, so steps is one less than input number
    r = r_min + (Dr/2.0)
    for i in range(r_steps-1):
        mu = mu_min + (Dmu/2.0)
        for j in range(mu_steps-1):
            nu = nu_min + (Dnu/2.0)
            for k in range(nu_steps-1):
                x,y,z = co.GC2xyz(mu, nu, r, params.wedge)
                rho = hernquist_profile(x,y,z,params.q,params.r0)
                if rho > max_rho:
                    max_rho = rho
                nu = nu + Dnu
                count = count + 1
                if count % (total_steps/10) == 0:
                    print "Progress: ", progress, "percent searched"
                    progress = progress + 10
            mu = mu + Dmu
        r = r + Dr
    return max_rho

def generate_perturbation(num_stars, params, parameters, batch=1000, fail_quit=100,
                          fileout="Pertgen82.txt", detection=1, convolve=1, append=0):
    """Density of perturbation added to background, as a function of position"""
    # Initialize file
    if append==0:
        out = open(fileout, 'w')
        out.write("# "+str(num_stars)+" stars, l,b,r \n")
        out.close()
    # Get integral
    tot_prob = get_max_prob_2(params, parameters)
    N_out, fails = 0, 0
    g_min, g_max = params.stripe[2][0], params.stripe[2][1]
    while N_out < num_stars:
        # Generate Points
        mu, nu, r = get_stripe_points(params.mu_lim, params.nu_lim, params.r_lim, batch)
        # Test points for inclusion
        x,y,z = co.GC2xyz(mu, nu, r, params.wedge)
        holder = []
        for i in range(len(mu)):
            rho = perturb_density(x[i],y[i],z[i], parameters)
            #print (rho / tot_prob), rho, tot_prob
            if (rho / tot_prob) > np.random.uniform():
                l,b,r1 = co.xyz2lbr (x[i], y[i], z[i])
                # Convolve
                if convolve==1:
                    r1 = star_convolution(r1)
                # Detection
                if detection==1:
                    m_g = co.getg(r1)
                    if np.random.uniform() > (sigmoid_error(m_g)):  continue
                if (co.getg(r1) < g_min) or (co.getg(r1) > g_max):  continue
                # Add to keepers
                holder.append([round(l,6),round(b,6),round(r1,6)])
                N_out = N_out + 1
        # Failure code
        if len(holder) == 0:  fails = fails + 1
        if fails >= fail_quit:  break
        # Remove possible excess stars
        if N_out > num_stars:
            slice = -1*(N_out-num_stars)
            holder = holder[:slice]
            print "#---Sliced {0} stars to make quota".format(str(-1*slice))
            N_out = N_out + slice #Slice is negative
        # Add to dataset
        if len(holder) != 0:
            if fi.append_data(sc.array(holder), fileout, delimiter=" ") == 1:
                print "#---Perturbation Progress: {0} stars of {1} total stars generated".format(N_out, num_stars)
    if fails >= fail_quit:  
        print "!!! Perturbation generation FAILED due to overstepping empty batch limit: {0}".format(fail_quit)
    else:
        print "#---Perturbation generation succeeded, written as {0}, with {1} empty batches".format(fileout, fails)
    return fileout

def get_max_prob_2(params, parameters):  #should be combined with get_max_prob()
    mu_min, mu_max, mu_steps = params.mu_lim
    nu_min, nu_max, nu_steps = params.nu_lim
    r_min, r_max, r_steps = params.r_lim
    total_steps = ((mu_steps-1)*(nu_steps-1)*(r_steps-1))
    Dmu = (mu_max-mu_min) / mu_steps
    Dnu = (nu_max-nu_min) / nu_steps
    Dr = (r_max-r_min) / r_steps
    print "Finding maximum probability of perturbation function,", total_steps, "steps"
    max_rho, count, progress = 0.0, 0, 0
    # Get value at center of each bin, so steps is one less than input number
    r = r_min + (Dr/2.0)
    for i in range(r_steps-1):
        mu = mu_min + (Dmu/2.0)
        for j in range(mu_steps-1):
            nu = nu_min + (Dnu/2.0)
            for k in range(nu_steps-1):
                x,y,z = co.GC2xyz(mu, nu, r, params.wedge)
                rho = perturb_density(x,y,z, parameters)  #Only change!!!
                if rho > max_rho:
                    max_rho = rho
                nu = nu + Dnu
                count = count + 1
                if count % (total_steps/10) == 0:
                    print "Progress: ", progress, "percent searched"
                    progress = progress + 10
            mu = mu + Dmu
        r = r + Dr
    return max_rho

def perturb_density(x,y,z, parameters):
    """Density due to added galaxy-centered error profile"""
    r = sc.sqrt(x*x + y*y + z*z)
    g = co.getg(r)
    exponent = (g*parameters[1]) + parameters[2]
    return (parameters[0] + sc.exp(exponent))
    
def sigmoid_error(x, modulus=None):
    """Application of detection efficiency"""    
    s = [0.9402, 1.6171, 23.5877]
    if modulus != None:
        s[2] = s[2] - modulus
    detection_efficiency = s[0] / (np.exp(s[1]*(x - s[2])) + 1.)
    return detection_efficiency

def star_convolution(r, mu=4.2, sigma=0.6):
    """ Convolve Stars based on turnoff distribution """
    m_g = co.getg(r)
    m_g = np.random.normal(m_g, 0.6)
    return co.getr(m_g)
    
def build_stripe(params, filename, num_stars, perturb_weight=None,
                 perturb_params=(0.0, 0.79, -19.9), con=1, det=1, app=1):
    """ Build a test SDSS data wedge, based on the parameters in the ParamSet; including
        background, perturbation, and streams """
    params.print_params()
    seed = int(time.time())
    np.random.seed(seed)
    print "Random Seed: {0}".format(seed)
    # If the perturbation is non-zero, assign the correct percentage of stars to 
    if perturb_weight != None:
        perturb_stars = int(num_stars*perturb_weight)
        num_stars = num_stars - perturb_stars
        print "# - {0} Stars in Perturbation".format(perturb_stars)
    # get weights
    denom = 1.0
    for i in range(len(params.streams)):
        denom = denom + np.exp(params.stream_weight[i])
    back_stars = int((1.0 / denom)*num_stars)
    print "# - {0} Stars in Background".format(back_stars)
    stream_stars = []
    for i in range(len(params.streams)):
        new_stars = int((np.exp(params.stream_weight[i]) / denom)*num_stars)
        stream_stars.append(new_stars)
        print "# - {0} Stars in Stream {1}".format(new_stars, i) 
    # app is zero here in order to generate file
    x = generate_background(back_stars, params, batch=10000, fail_quit=100,
                        fileout=filename, detection=det, convolve=con, append=0)
    print "# --- Background Generation Complete"
    if len(stream_stars) > 0:
        for i in range(len(stream_stars)):
            y = stream_into_stripe(params, i, stream_stars[i], batch=10000, fileout=filename,
                        detection=det, convolve=con, append=app)
            print "# --- Stream {0} Generation Complete".format(i)
    if perturb_weight != None:
        z = generate_perturbation(perturb_stars, params, perturb_params, batch=1000,
                        fail_quit=100, fileout=filename, detection=det, convolve=con, append=app)
        print "# --- Perturbation Generation Complete"
    return 1

'''Need output in l, b, r, with number of stars at the top '''

if __name__ == "__main__":
    params = ParamSet()
    build_stripe(params, "detcon_test_82_1.txt", 100000, perturb_weight=0.05,
                 perturb_params=(0.0, 0.79, -19.9), con=1, det=1, app=1)
    params.print_params()
    print "\n ### Done"

""" ---------------------  Not used --------------------------------- """
def rough3d_integrator(params):  #q, r0, mu_lim, nu_lim, r_lim):
    """ lims are (min, max, N_steps) """
    mu_min, mu_max, mu_steps = params.mu_lim
    nu_min, nu_max, nu_steps = params.nu_lim
    r_min, r_max, r_steps = params.r_lim
    total_steps = ((mu_steps-1)*(nu_steps-1)*(r_steps-1))
    Dmu = (mu_max-mu_min) / mu_steps
    Dnu = (nu_max-nu_min) / nu_steps
    Dr = (r_max-r_min) / r_steps
    print "Integrating function,", total_steps, "steps"
    integral, count, progress = 0.0, 0, 0
    # Integrate from value at center of each bin, so steps is one less than input number
    r = r_min + (Dr/2.0)
    for i in range(r_steps-1):
        mu = mu_min + (Dmu/2.0)
        for j in range(mu_steps-1):
            nu = nu_min + (Dnu/2.0)
            for k in range(nu_steps-1):
                #print mu, nu, r
                x,y,z = co.GC2xyz(mu, nu, r, params.wedge)
                rho = hernquist_profile(x,y,z,params.q,params.r0)
                dV = r*r*sc.cos(nu)*Dmu*Dnu*Dr
                integral = integral + rho*dV
                nu = nu + Dnu
                count = count + 1
                if count % (total_steps/20) == 0:
                    print "Progress: ", progress, "percent complete"
                    progress = progress + 5
            mu = mu + Dmu
        r = r + Dr
    return integral
