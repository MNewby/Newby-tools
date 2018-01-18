import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func

test_bastard

'''python script for quickly analyzing data.
Matthew Newby, November 15,2010'''

x = sc.arange(0.0, 100.0, 0.1)
#y = func.test_bastard(x, params)
#y = func.plank_curve(x, [ 0.001, 1.0, 0.0])
#y = func.synaptic_response(x, [10.0])
#y = func.sigmoid_function(x, [1.0, 10.0, 0.0])
y =  (1.0/(1.0 + np.exp( (x-10.0)/(-5.0) ))) * (1.0/(1.0 + np.exp( (x-15.0)/(40.0) )))
plt.figure()
plt.plot(x,y)
plt.show()
plt.close('all')
