import math as ma
import numpy as np
import scipy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt

"""python script 
Matthew Newby, RPI June 6, 2011"""

a = 1.08574

def usual_mag (f, m0):
    return (m0 - 2.5*sc.log(f))

def var_mag (f, sigma):
    return (a*a*sigma*sigma)/(f*f)
    
def luptitude(f, b, mu0):
    holder = a*( sc.arcsinh(f/(2.0*b)) )
    return (mu0 - holder)
    
def var_lup(f, sigma, b):
    return (a*a*sigma*sigma)/( (4.0*b*b) + (f*f))
    
if __name__ == "__main__":
    sigma = 10.0
    b = 0.53 #1.042*sigma
    f0 = 20.0
    m0 = 24.22 #2.5*sc.log(f0)
    mu0 = 24.91 #(m0 - 2.5*sc.log(b))
    f = sc.arange(0.0, 40.0, 1.0)
    y1 = usual_mag(f, m0)
    y2 = luptitude(f, b, mu0)
    y3 = var_mag(f, sigma)
    y4 = var_lup(f, sigma, b)
    plt.figure()
    #plt.plot(f,y1)
    #plt.plot(f,y2)
    print "Magnitude =", usual_mag(3.0, m0)
    print "Luptitiude = ", luptitude(3.0, b, mu0)
    plt.plot(f,y3)
    plt.plot(f,y4)
    plt.show()