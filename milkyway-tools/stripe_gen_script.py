import math as ma
import numpy as np
import scipy as sc
import time
import test_wedge_generator as twg

'''python script for quickly building several test data wedges.
Matthew Newby, November 1, 2011'''

# Set stream generation parameters
""" #Nate's
background = [0.0, 1.0, 0.670, 13.50, 1.0]
stream1 = [-1.828, 31.361, 29.228, 1.445, 3.186, 2.854]
stream2 = [-2.182, 1.0, 25.0, 1.0, -1.660, 1.50]
"""

# Mine, simpler
background = [0.0, 1.0, 0.70, 19.50, 1.0]
stream1 = [-2.0, 1.0, 14.5, 1.5, 3.0, 2.0]
stream2 = [-2.5, 40.0, 23.0, 1.0, -1.660, 1.50]

single = twg.ParamSet()
single.background = background
single.streams = [stream1]

double = twg.ParamSet()
double.background = background
double.streams = [stream1, stream2]
double.update_refs()

perturbation = [None] #, 0.05, 0.10, 0.15, 0.20]
per_params = (0.0, 0.79, -19.9)

t1 = time.time()

"""
# Single Stream Generation
for i in range(len(perturbation)):
    t2 = time.time()
    filename = "single_t82_"+str(i)+"_new.txt"
    print "### --- Starting run {0} with {1} stars in perturbation".format(filename, perturbation[i])
    twg.build_stripe(single, filename, 110000, perturbation[i], per_params,
                 con=1, det=1, app=1)
    t3 = time.time()
    print "# - completed in {0} senconds".format(t3 - t2)

print " * single stream runs completed in {0} seconds".format(t3-t1)
"""

# Double Stream Generation
for i in range(len(perturbation)):
    t2 = time.time()
    filename = "double_t82_"+str(i)+"_new.txt"
    print "### --- Starting run {0} with {1} stars in perturbation".format(filename, perturbation[i])
    twg.build_stripe(double, filename, 110000, perturbation[i], per_params,
                 con=1, det=1, app=1)
    t4 = time.time()
    print "# - completed in {0} seconds".format(t4-t2)

#print " * double stream runs completed in {0} seconds".format(t4-t3)
print " * all runs completed in {0} seconds".format(t4-t1)
print " --- fin --- "
