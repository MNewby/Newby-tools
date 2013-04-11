import math as ma
import numpy as np
import scipy as sc
import astro_coordinates as coor
import sgr_law as sl

'''python script for figuring out what l,b values were used to get a Sgr Lambda, Beta
Matthew Newby, July 16, 2012'''

pres = 0.001  #search precision

center_l, center_b, r = 180.0, 0.0, coor.getr(20.75)  #search start

lam_out = [75.0, 85.0, 95.0, 105.0, 115.0]
#beta_out = [9.5, 10.3, 10.3, 10.7, 11.8]  # Secondary Peak
beta_out = [-1.3, -2.5, -1.3, -0.9, -1.4]  # Primary Peak

l_steps = sc.arange(-180.0, 190.0, 20.0)  # -100 to 100, increments of 20
b_steps = sc.arange(-90.0, 91.0, 10.0)  # -20 to -20, increments of 4

steps = len(l_steps)


def dist(l,b,lam, beta):
    x,y,z = coor.lbr2xyz(l,b,r)
    x0,y0,z0,lam0,beta0,r0 = sl.law_xyz2sgr_sun(x,y,z)
    one = (lam0 - lam)*(lam0 - lam)
    two = (beta0 - beta)*(beta0 - beta)
    return sc.sqrt(one + two)

verbose = 0

for i in range(len(lam_out)):
    d0, scale = 100.0, 1.0
    l0, b0 = 0.0, 0.0
    l_base, b_base = center_l, center_b
    lam, beta = lam_out[i], beta_out[i]
    x, changed = 0, 0
    while x<15:
        for dl in l_steps:
            l = l_base + dl*scale
            for db in b_steps:
                b = b_base + db*scale
                d = dist(l,b,lam,beta)
                if d < d0:  l0, b0, d0 = l, b, d; changed=1
        if verbose==1:
            print "Current l,b,d:  {0}, {1}, {2}".format(l0,b0,d0)
            x0,y0,z0 = coor.lbr2xyz(l0,b0,r)
            print sl.law_xyz2sgr_sun(x0,y0,z0)[3:]
        if d0 < pres:  print "# - Exiting by threshold"; break
        scale = scale/2.0
        l_base, b_base = l0, b0
        x=x+1
    print "# - For lamda, beta = ({0},{1}), l,b,d = {2}, {3}, {4}".format(lam,beta,l0,b0,d0)
    x0,y0,z0 = coor.lbr2xyz(l0,b0,r)
    print "# - {0}".format(sl.law_xyz2sgr_sun(x0,y0,z0)[3:])