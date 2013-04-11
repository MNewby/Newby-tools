import math as ma
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import astro_coordinates as co

rad = ma.pi/180.0
deg = 180.0/ma.pi

G = 6.67E-11
M = 1.99E30

r = sc.arange(0.1, 70.1, 0.1)
rm = r*1.496E11  #AU to meters
vm = sc.sqrt(G*M/rm)  #velocity in m/s
v = vm/1000.0  # velocity in km/s

#planets
np = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Eris", "Ceres"]
rp = [0.467, 0.728, 1.000, 1.666, 5.204, 9.582, 19.229, 30.104, 39.264, 68.01, 2.765]  # radius, AU
vp = [47.87, 35.02, 29.78, 24.08, 13.07, 9.69,  6.81,   5.43,   4.70,   3.43,  17.88]  # radial v, km/s
sp = sc.array([0.383, 0.950, 1.000, 0.533, 11.2,  9.45,  4.01,   3.88,   0.18,   0.188, 0.076])  # radius in Earth rad
cp = ["grey", "green", "blue", "red", "orange", "yellow", "cyan", "cyan", "black", "black", "grey"]

""" Plot 1 - Solar System Rotation Curve
plt.figure(1)
plt.plot(r,v, "k--")
plt.scatter(rp, vp, c=cp, marker="o", s=75*(sc.log(sp)+2.0))
for i in range(len(np)):
    if     np[i]=="Eris":  plt.text(rp[i]-1.0, vp[i]+0.5, np[i], fontsize=10)
    else:  plt.text(rp[i]+1.0, vp[i]+0.5, np[i], fontsize=10)
plt.fill([30.0, 65.0, 65.0, 30.0],[45.0, 45.0, 30.5, 30.5], 'b', alpha=0.25)
plt.text(32.0, 42.0, "Solar System Rotation Curve", fontsize=14)
plt.text(47.5, 32.0, "Point size corresponds to log$_{10}$ planet radii\n (Credit: Matthew Newby, Milkyway@home)", 
    horizontalalignment="center", fontsize=10)
#\n (Matthew Newby, Milkyway@home)")
plt.text(40.0, 37.0, r"$v = \sqrt{GM/r}$", fontsize=20)
plt.xlim(0.0, 70.0)
plt.ylim(0.0, 50.0)
plt.xlabel(r"Radius, $r$ (AU)")
plt.ylabel(r"Orbital Velocity, $v$ (km/s)")
"""

""" Plot number 2 - Galactic rotation curves """
plt.figure(2)
hh = 1.0
aa = 1.0
vv = 3.0
dd = 13.0
x = sc.arange(0.01, 15.01, 0.01)
disk = (hh*(1.0 - sc.exp(-x/hh))/x)#*(1.0 / (1.0 + sc.exp(2.0*(aa-x))))
plummer = ( (x*x) / ((x*x + aa*aa)**(3.0/2.0)) )
halo = vv*vv*((2.0*x) / (x*x + dd*dd))
y2 = sc.sqrt(disk)
y3 = 0.8*sc.sqrt(plummer)
y4 = sc.sqrt(halo)
y1 = y2 + y3
y0 = y2 + y3 + y4
plt.plot(x,y1, linewidth=2.0)
plt.plot(x,y0*0.84, linewidth=2.0)
plt.plot(x,y2, "k--")
plt.plot(x,y3, "k-.")
plt.plot(x,y4, "k:")
plt.text(11.7, 1.16, "Observed", fontweight="bold")
plt.text(11.7, 0.55, "Expected", fontweight="bold")
plt.text(12.3, 0.85, "Dark Matter")
plt.text(12.3, 0.31, "Disk (light)")
plt.text(12.3, 0.14, "Bulge (light)")
plt.text(10.0, 1.5, "Galaxy Rotation Curve", fontsize=16, horizontalalignment="center")
plt.text(10.0, 1.44, "Credit:  Matthew Newby, Milkyway@home", fontsize=10, horizontalalignment="center")
plt.xlim(0.0, 15.0)
plt.ylim(0.0, 1.6)
plt.xlabel(r"Radius, $r$ (kpc)")
plt.ylabel(r"Orbital Velocity ($arb.$)")
plt.show() 
