import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func
import time

'''python script for quickly creating test data sets.
Matthew Newby, May 26, 2011'''

out_file = "test_poly_out.txt"
function = func.Nth_order_polynomial
parameters = [11.0, -1.0, 0.5]

x = sc.arange(-4.0, 4.0, 0.1)
y = function(x, parameters)
sig = 1.0  #1 sigma value for errors, takes:  None or float 

out = sc.ones((len(x), 3))
out[:,0] = x
if (sig == None):  out[:,1] = y
else:
    np.random.seed( int(time.time()) )
    for i in range(len(y)):
        out[i,1] = np.random.normal(y[i], sig)
        out[i,2] = float(sig)

if fi.write_data(out, out_file, header="#-x, y, sig"):
    print "#- File output succeeded"
    
plt.figure(1)
plt.scatter(out[:,0],out[:,1],2,'k','o')
plt.show()
plt.close('all')
print "# - Done"
    