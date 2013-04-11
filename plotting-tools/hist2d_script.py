import math as ma
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import astro_coordinates as co
import files as fi


import plot_pack as pp

#data1 = fi.read_data("NGC_5272.csv", ",")
data2 = fi.read_data("NGC_6205.csv", ",")

x = data2[:,2]-data2[:,3]
y = data2[:,2]

gerr = 3.11E-4 + np.exp(0.79*data2[:,2] - 20.0)
rerr = -2.63E-5 + np.exp(0.80*data2[:,3] - 19.8)
xerr = np.sqrt(gerr*gerr + rerr*rerr)
yerr = gerr

#Initializations
tester = pp.HistMaker(x, y, xsize=0.01, ysize=0.1, xarea=None, yarea=(15.0, 25.0) )
tester.cmap = 'color'             # Sets colormap
tester.varea = (0.0,None)         # Sets contrast range for pixels
tester.scale = None               # Special scaling, such as 'sqrt' or 'log'
tester.labels = [r'$(g-r)_0$',r'$g$']      # Labels for x,y axes
tester.xflip = 0   # Set=1 to flip respective axis
tester.yflip = 1   # Set=1 to flip respective axis

#tester.gblur(blur=1)  # Uncomment this to apply a Gaussian blur to the image

# Use this next line to bin 
#tester.H = pp.ErrBlurPoint(x, y, tester, con=30, xerror=xerr, yerror=yerr)

tester.plot()

"""                
import scipy.special as ss
sq2 = ma.sqrt(2.0)


def getb(sig, a, p):
    part1 = ss.erf( a / (sq2*sig) )
    part2 = 2.0*p + part1
    part3 = sq2*sig*ss.erfinv(part2)
    return part3

mu = 0.0
sig = 1.0
con = 120

p = 0.9973 / float(con)

xpoints = [0.0]
for i in range(con/2):  xpoints.append(getb(sig, xpoints[i], p) )
x = -1.0*sc.array(xpoints[1:])
xx = sc.append(sc.array(xpoints), x)

ypoints = [0.0]
for i in range(con/2):  ypoints.append(getb(sig, ypoints[i], p) )
y = -1.0*sc.array(ypoints[1:])
yy = sc.append(sc.array(ypoints), x)

newx = []
newy = []
for i in range(len(xx)):
    for j in range(len(yy)):
        newx.append(xx[i])
        newy.append(yy[j])
surface, xedge, yedge = np.histogram2d(newx, newy, 21)
X = xedge
Y = yedge
#if len(X) > surface.shape[1]:  X = X[:-1]  # Needed in case of rounding error
#if len(Y) > surface.shape[0]:  Y = Y[:-1]  # Needed in case of rounding error

print surface.shape, len(X), len(Y)

plt.figure(1)
#plt.hist(xx, 21, align='mid')
plt.scatter(newx, newy, s=1, facecolor="black")
#plt.imshow(np.sqrt(surface), origin='lower')
plt.contour(X[:-1], Y[:-1], surface)
plt.show()
"""



"""
plt.figure(1)
#sp1 = plt.subplot(121)
#plt.scatter(data1[:,2]-data1[:,3], data1[:,2], s=1)
#sp2 = plt.subplot(122, sharex=sp1, sharey=sp1)
plt.scatter(data2[:,2]-data2[:,3], data2[:,2], s=1)
plt.ylim(25.0, 15.0)
plt.xlim(0.0, 1.0)
plt.show()
"""
