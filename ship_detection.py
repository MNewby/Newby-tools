import math as ma
import numpy as np
import scipy as sc
from scipy import integrate
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import astro_coordinates as co

rad = ma.pi/180.0
deg = 180.0/ma.pi

""" Script for determining the brightness magnitude, angular resolution and proper 
    motion of a theoretical spaceship of face-on surface area 'A', distance 'd', 
    temperature 'T', and wavelength range 'l_min' to 'l_max'.  
    Telescope radius 'tele_r'
    
    <Maybe add options for thrust and other heat generation>
    
    -Matthew Newby, Jan 03, 2013
    """

h = 6.62606957E-34  # Plank's Constant [J*s]
k = 1.3806488E-23   # Boltzmann's Constant [J/K]
c = 2.99792458E8    # Speed of Light [m/s]
#tele_A = 2.0*ma.pi*tele_r*tele_r  #NOT NEEDED

def get_mag(A, d, T, l_min, l_max, tele_r):
    sigma = 5.67E-8
    D = d*4.85E-6  # AU to parsecs  UNITS???
    total_power_out = A*sigma*T*T*T*T  # Stephan-Boltzmann law [Watts]
    # multiplies by fraction of spectrum actually detected
    power_out = get_frac(T, l_min, l_max)*total_power_out
    print "{0} Watts detected out of {1} Watts emitted".format(power_out, total_power_out)
    flux = power_out / (4.0*ma.pi*D*D)  # FIX THIS! UNITS?
    flux = flux / 10E-26  #to Janskys - UNITS?
    print flux
    # See table 1.2 in S&G!
    mag = -2.5*np.log10(flux)+8.90-0.0  #ZP
    return mag

def get_frac(T, l_min, l_max):
    nu_min, nu_max = c/l_max, c/l_min
    num, err = integrate.quad(planck_law, nu_min, nu_max, T)
    frac = (num*15.0)/(ma.pi**4)
    return frac
    
def planck_law(nu, T):
    pre = (2.0*h*nu*nu*nu) / (c*c)
    ex = ma.exp( (h*nu)/(k*T))
    return pre*(1.0/(ex-1.0))


if __name__ == "__main__":
    A = 10.0        #(in m^2)
    d = 1.0         #(in AU)
    T = 300.0       #(in Kelvin)
    l_min = 507.0   #(in nm)
    l_max = 595.0   #(in nm)
    tele_r = 0.005  #(in m)
    print get_mag(A, d, T, l_min, l_max, tele_r)
