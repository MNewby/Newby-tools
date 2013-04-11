import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import files as f

#No longer used - has been encapsulated and included in gctools_enc.py

#d = raw_input('GC distance:')
#distance = eval(d)
#starfile = raw_input('Star file name:')
original_file = "NGC_7089_g-r.dat"
convolved_file = "convolved_7089_23.2kpc_g-r.dat"

left_limit = 0.0
right_limit = 0.3

#unconvolved analysis
y = f.read_data(original_file)
l, w = y.shape

o_right = 0
o_left = 0
o_middle = 0
for i in range(l):
    if (y[i,1] < left_limit):
        o_left = o_left + 1
    elif (y[i,1] > right_limit):
        o_right = o_right + 1
    else:
        o_middle = o_middle + 1
print "Before Convolution:"
print "Left:", o_left, "in region:", o_middle, "Right:", o_right


#Convolved Analysis
z = f.read_data(convolved_file)
l, w = z.shape

c_right = 0
c_left = 0
c_middle = 0
for i in range(l):
    if (z[i,1] < left_limit):
        c_left = c_left + 1
    elif (z[i,1] > right_limit):
        c_right = c_right + 1
    else:
        c_middle = c_middle + 1
print "After convolution:"
print "Left:", c_left, "in region:", c_middle, "Right:", c_right

#Finding differences
delta_l = o_left - c_left      #These are fluxes, so they are reversed
delta_m = c_middle - o_middle  #This is a standard "final-minus-initial"
delta_r = o_right - c_right    #These are fluxes, so they are reversed

print "Flux into f-turnoff region, from left:", delta_l
print "Flux into f-turnoff region, from right:", delta_r
print "Total change in f-turnoff stars:", delta_m
