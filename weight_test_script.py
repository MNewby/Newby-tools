import math as ma
import numpy as np
import scipy as sc
import files as fi

'''python script for quickly analyzing distribution function weights.
Matthew Newby, November 15, 2011'''

N_stars = 1000
s1_w = -2.0
s2_w = -5.0
pp = 0.10

p_w = np.log((1.0 + np.exp(s1_w) + np.exp(s2_w)) / ( (1.0/pp) - 1.0))

denom = 1.0 + np.exp(s1_w) + np.exp(s2_w) + np.exp(p_w)

N_back = int((1.0 / denom)*N_stars)
N_s1 = int((np.exp(s1_w)/denom)*N_stars)
N_s2 = int((np.exp(s2_w)/denom)*N_stars)
N_p = int((np.exp(p_w)/denom)*N_stars)

print N_back, N_s1, N_s2, N_p, (N_back+N_s1+N_s2+N_p)