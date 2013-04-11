import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import astro_coordinates as co

rad = ma.pi/180.0
deg = 180.0/ma.pi


def NFW(r, a=10.0):
    """ Gives acceleration in an NFW halo"""
    term1 = ( a/(r*r) ) * np.log(1.0 + (r/a) )
    term2 = (r * (1.0 + (r/a) ) )
    # Prefix:  -4.0*ma.pi*G*rho_0
    return (term1 - (1.0/term2))*a*a
    
def log_halo(x,y,z, q=0.9, d=13.0, Vh=220.0):
    """ Gives acceleration in an log-halo, values from Law+2005 Oblate halo fit"""
    top = 2.0*np.sqrt(x*x + y*y) + 2.0*(z/(q*q))
    bottom = (x*x + y*y) + ((z*z)/(q*q)) + (d*d)
    return Vh*Vh*top/bottom
    
     
caustic_x = fi.read_data("constzyrangex2.txt")
caustic_y = fi.read_data("constxzrangey2.txt")
caustic_z = fi.read_data("constxyrangez2.txt")

radius = np.arange(1.0, 256.01, 0.1)
acc = NFW(radius)

logx, logz = np.arange(1.0, 256.01, 0.1), np.arange(1.0, 256.01, 0.1)
flat = np.zeros(len(logx), float)
alogx = log_halo(logx, flat, flat)
alogz = log_halo(flat, flat, logz)


norm1 = 1.0/900.0 # np.ma.max(acc)
norm2 = 1.0 #np.ma.max(caustic_x)
norm3 = 1.0 #np.ma.max(caustic_y)
norm4 = 1.0 #np.ma.max(caustic_z)
norm5 = 4.0 #np.ma.max(caustic_z)
norm6 = 4.0 #np.ma.max(caustic_z)

plt.figure(1)
plt.plot(radius, acc/norm1, label="NFW")
plt.plot(caustic_x[:,0], caustic_x[:,2]/norm2, label="X")
plt.plot(caustic_y[:,0], caustic_y[:,2]/norm3, label="Y")
plt.plot(caustic_z[:,0], caustic_z[:,2]/norm4, label="Z")
plt.plot(logx, alogx/norm5, label="logx")
plt.plot(logz, alogz/norm6, label="logz")
plt.legend(loc='upper right')
plt.xlim(0.0, 260.0)
plt.ylim(0.0, 1200.0)
plt.xlabel("Distance (kpc)")
plt.ylabel(r"acceleration (kpc/Gyr${}^2$)")
plt.show()
