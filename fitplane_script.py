import math as ma
import numpy as np
import scipy as sc
import functions as func
import files as fi
import fit as fit
import plane_fit as pf
import astro_coordinates as co

#data = fi.read_data("../sgrnorth_paper/Newberg2007_Sgr.txt", ",")
#data = fi.read_data("../sgrnorth_paper/Koposov2012_Sgr.txt", ",")
data = fi.read_data("../sgrnorth_paper/Koposov2012_2nd.txt", ",")
# In fit parameter format:
#data = fi.read_data("../sgrnorth_paper/sgr_params_south.txt", ",")
#data = fi.read_data("../sgrnorth_paper/sgr_params_north.txt", ",")

fitter=fit.ToFit(data[:,0], data[:,1], "1.0")
"""
# gots to be tricksy, yes...  plane fitting requires dark magiks...
#wedges = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
#wedges = [79, 82, 86]
# errors from sgr_north_utils.py
#errors = [9.31738638849, 9.1020519394, 19.7856077035, 20.8318952656, 9.39314333747,
#          26.9966491404, 40.6343350862, 15.3631542798, 11.270906999, 10.3242787823,
#          10.290985838, 8.75140991672, 14.6281676367, 16.0532641683, 6.69586402716]
#errors = [10.3012569286, 9.26151305927, 10.4476312495]  # Errors are wrong!!
# These are right?
# 7.57640162772, 7.67810980299, 11.5363484879, 25.6937588232, 8.75918164562,
# 37.6255221984, 54.5484852169, 13.1237001504, 12.2803690866, 11.9300143607,
# 12.2150368604, 13.3084651939, 14.9866923184, 15.7846605126, 27.9365850481, 
# 11.2583497941, 10.1111224829, 9.12819009006
to_x = sc.zeros((len(data[:,0]), 3) )  # RAISE THIS TO 4? WHEN ERRORS ARE USED
for i in range(len(data[:,0])):
    u, v, w, mu, r, theta, phi, sigma = 0.0, 0.0, 0.0, data[i,4], data[i,5], \
        data[i,6], data[i,7], data[i,8]
    to_x[i,0], to_x[i,1], to_x[i,2] = co.stream2xyz(u, v, w, mu, r, theta, phi, wedges[i])
    #to_x[i,3] = errors[i]

"""

to_x = sc.zeros((len(data[:,0]), 3) )
for i in range(len(data[:,0])):
    x,y,z = co.lbr2xyz(data[i,0], data[i,1], co.getr(data[i,2], 0.6))  #<--- 0.6 for Koposov 'i's
    to_x[i,0], to_x[i,1], to_x[i,2] = x,y,z

print to_x


params = pf.plane_OLS(to_x[:,0], to_x[:,1], to_x[:,2])
print "#- Distance array:", pf.plane_dist(to_x[:,0], to_x[:,1], to_x[:,2], params)

fitter.x = to_x
fitter.y = sc.zeros(len(data[:,0]) )
fitter.function=func.plane_fit  #func.R_squared_plane
fitter.update_params([3.0*ma.pi/2.0, ma.pi/2.0, 1.0]) #(params)
fitter.range = [(0, 0.0, 2.0*ma.pi, ma.pi/10.0), (1, 0.0, ma.pi, ma.pi/20.0), (2, -20.0, 20.0, 1.0)]
print fitter.params, fitter.RR

fitter.step = [0.1, 0.1, 0.1]
path1 = fit.MCMC(fitter, iters=10000, annealing=0, verbose=0)

fitter.update_params([3.0*ma.pi/2.0, ma.pi/2.0, 1.0])
fitter.step = [0.75, 0.75, 0.75]
print "# - Steps: ", fitter.step
path2 = fit.gradient_descent(fitter, its=10000, line_search=0)
errors = fit.get_Hessian_errors(fitter)

#fitter.plot()
yi, xi = 1,0
surface = fit.get_surface(fitter, yi, xi, plot=1, path=(path1[:,yi], path1[:,xi]))
surface = fit.get_surface(fitter, yi, xi, plot=1, path=(path2[:,yi], path2[:,xi]))
#print fitter.params, fitter.RR

t, p, d = fitter.params
dt, dp, dd = errors
params = [sc.cos(t)*sc.sin(p), sc.sin(t)*sc.sin(p), sc.cos(p), d]
errs = [(-sc.sin(t)*sc.sin(p)*dt) + (sc.cos(t)*sc.cos(p)*dp),
        (sc.cos(t)*sc.sin(p)*dt) + (sc.sin(t)*sc.cos(p)*dp),
        (-sc.sin(p)*dp), dd]
bottom = sc.sqrt(params[0]*params[0] + params[1]*params[1] + params[2]*params[2])
print bottom
print params
print errs

# Correlation:

# x = distance along plane
x0 = to_x[0,:]
XX = sc.sqrt((to_x[:,0]-x0[0])**2 + (to_x[:,1]-x0[1])**2 + (to_x[:,2]-x0[2])**2)
# y = dustance from plane
YY = pf.plane_dist(to_x[:,0], to_x[:,1], to_x[:,2], params)
#YY = func.plane_fit(fitter.x, fitter.params)
print "Correlation r: {0}".format(fit.get_correlation(XX,YY))

"""
# Plot of distance from plane
import matplotlib
import matplotlib.pyplot as plt

plt.figure()
#YY = func.plane_fit(fitter.x, fitter.params)
#XX = fitter.x[:,0] #wedges
plt.scatter(XX,YY)
plt.show()

#print sc.sum(YY)
#print sc.std(YY)
#print sc.sum(YY)/sc.std(YY)
"""
