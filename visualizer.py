import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import test_wedge_generator as twg
import astro_coordinates as co

'''python script for visualizing test data stripe.
Matthew Newby, July 7, 2011'''

"""
#Builds and displays a comparison of data density compared with herquist model
#data = fi.read_data("Backgen82.txt")
data = fi.read_data("Pertgen82.txt")
#x, y, z = co.lbr2xyz(data[:,0], data[:,1], data[:,2])
#d = sc.sqrt(x*x + y*y + z*z)
r = data[:,2]
size = 1.0
max, min = np.ma.max(r), np.ma.min(r)
bins = int((max - min + 1.0) / size)
rho_raw = sc.zeros(bins, int)
for i in range(len(r)):
    index = int( (r[i]-min) / size)
    rho_raw[index] = rho_raw[index] + 1
# get volume at each distance
#r_min, r_max = co.getr(16.0), co.getr(22.5)
mu_min, mu_max = 310.0, 419.0
nu_min, nu_max = (-1.25*ma.pi/180.0), (1.25*ma.pi/180.0)
mu_bit = mu_max - mu_min
nu_bit = sc.sin(nu_max) - sc.sin(nu_min)
V = sc.zeros(bins, float)
for i in range(bins):
    r_bit = (min + (i+1)*size)**3 - (min + i*size)**3
    V[i] = r_bit*(1./3.)*mu_bit*nu_bit
# make the histogram of rho(r)
rho_hist = sc.zeros((bins,2), float)
for i in range(bins):
    rho_hist[i,0] = min + (i*size)
    rho_hist[i,1] = rho_raw[i] / V[i]
#get gal-centered herquist distribution for r points in wedge
#mu_avg = ((mu_max - mu_min) / 2.0) + mu_min
#nu_avg = 0.0
#sol_r = sc.arange(min, max, (size/2.0))
#x,y,z = co.GC2xyz(mu_avg, nu_avg, sol_r, 82)
#gal_rho = twg.hernquist_profile(x,y,z,0.455,19.5)*1000000.0
#print gal_rho
#get gal-centered perturbation distribution for r points in wedge
mu_avg = ((mu_max - mu_min) / 2.0) + mu_min
nu_avg = 0.0
sol_r = sc.arange(min, max, (size/2.0))
x,y,z = co.GC2xyz(mu_avg, nu_avg, sol_r, 82)
gal_rho = twg.perturb_density(x,y,z, (0.0, 0.79, -19.9) )*0.78
print gal_rho
#Plot it all
fig = plt.figure()
plt.bar(rho_hist[:,0], rho_hist[:,1], size)
plt.plot(sol_r, gal_rho, c='r')
plt.xlabel(r'Sun-centered r (kpc)')
plt.ylabel(r'density')
plt.show()
"""


# Plotting Code 
fig = plt.figure()
ax = Axes3D(fig)

#ax.plot([-5.0,5.0],[0.0, 0.0],[0.0, 0.0], c='yellow')
#data_b = fi.read_data("Backgen82.txt")
#xb, yb, zb = co.lbr2xyz(data_b[:,0], data_b[:,1], data_b[:,2])
#ax.scatter(xb,yb,zb, facecolor='k', edgecolor='k')

data_p = fi.read_data("full_test_82.txt")
xp, yp, zp = co.lbr2xyz(data_p[:,0], data_p[:,1], data_p[:,2])
ax.scatter(xp,yp,zp, c='r')

#data_s = fi.read_data("stream_test.txt")
#xs, ys, zs = co.lbr2xyz(data_s[:,0], data_s[:,1], data_s[:,2])
#u,v,w = twg.generate_stream(1000, -2.2, 1.9, 2.0)
#xs,ys,zs = [],[],[]
#for i in range(len(u)):
#    x_h,y_h,z_h = co.stream2xyz(u[i],v[i],w[i],360.0, 21.0, 0.0, 0.0, 82)
#    xs.append(x_h), ys.append(y_h), zs.append(z_h)
#ax.scatter(xs,ys,zs, facecolor='b', edgecolor='b')
#mu, nu, r1 = sc.array([360.0, 360.0]),sc.array([-1.25, -1.25]),sc.array([10.0, 30.0])
#x,y,z = co.GC2xyz(mu, nu, r1, 82)
#ax.plot(x,y,z, c='r')
#mu2, nu2, r2 = sc.array([360.0, 360.0]),sc.array([1.25, 1.25]),sc.array([10.0, 30.0])
#x2,y2,z2 = co.GC2xyz(mu2, nu2, r2, 82)
#ax.plot(x2,y2,z2, c='r')
#xu1, yu1, zu1 = co.stream2xyz(10.0, 0.0, 0.0, 360.0, 21.0, 0.0, 0.0, 82)
#xu2, yu2, zu2 = co.stream2xyz(-10.0, 0.0, 0.0, 360.0, 21.0, 0.0, 0.0, 82)
#ax.plot([xu1,xu2],[yu1,yu2],[zu1,zu2],c='g')
#cset = ax.contour(x, y, z, zdir='z', offset=-10)
#cset = ax.contour(X, Y, Z, zdir='x', offset=-5.0)
#cset = ax.contour(X, Y, Z, zdir='y', offset=10)
ax.set_xlabel('x')
#ax.set_xlim3d(-40.0, 20.0)
ax.set_ylabel('y')
#ax.set_ylim3d(-5.0, 30.0)
ax.set_zlabel('z')
#ax.set_zlim3d(-40.0, 0.0)
plt.show()
