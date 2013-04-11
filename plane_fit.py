import math as ma
import numpy as np
from numpy import linalg as LA
import scipy as sc

def plane_OLS(x,y,z):
    """ Solves for the best-fit plane to a set of x,y,z data using ordinary
        least-squares.  Equation is of form z = Ax + By + C.
        Output is normalized a,b,c,d of a plane of the form ax+by+cz+d=0"""
    sum_x, sum_y, sum_z = sc.sum(x), sc.sum(y), sc.sum(z)
    sum_xx, sum_yy, sum_n = sc.sum(x*x), sc.sum(y*y), len(x)
    sum_xy, sum_xz, sum_yz = sc.sum(x*y), sc.sum(x*z), sc.sum(y*z)
    #  Will solve Ax=B, x is solution parameters (called p below)
    A = sc.matrix([ [sum_xx, sum_xy, sum_x],
                    [sum_xy, sum_yy, sum_y],
                    [sum_x,  sum_y,  sum_n] ])
    B = sc.matrix([sum_xz, sum_yz, sum_z])
    #print LA.inv(A)
    #print B.T
    p = LA.inv(A)*B.T
    params = [-float(p[0]), -float(p[1]), 1.0, -float(p[2])]  #c=1.0 by default, minuses due to OLS definition
    bottom = sc.sqrt(params[0]*params[0] + params[1]*params[1] + params[2]*params[2])
    for i in range(len(params)):  params[i] = params[i]/bottom
    print "# - Normalized best-fit plane parameters: {0}".format(params)
    return params

def plane_dist(x,y,z, params):
    a,b,c,d = params
    return (a*x + b*y + c*z + d)
    
if __name__ == "__main__":
    x = sc.array([1.0, 5.0, 0.0])
    y = sc.array([2.0, 2.0, 2.0])
    z = sc.array([-11.0, -51.0, -1.0])
    plane_OLS(x,y,z)