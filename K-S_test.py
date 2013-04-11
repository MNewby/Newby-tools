import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func
import time as ti

'''python script for testing the statistical strength of a distribution against
a null distribution.  Uses Kolmogorov-Smirnov test.
Matthew Newby, November 2, 2011'''

def get_CDF(unsorted_data):
    """ Generates the cumulative distribution function for a 1D set of data"""
    N = len(unsorted_data)
    cdf_x = np.sort(unsorted_data)  #The sorted data is also the x-axis
    cdf_y = sc.arange((1.0/N), (1.0+(1.0/N)), (1.0/N))
    if len(cdf_x) != len(cdf_y):  print "!!! error buidling CDF !!!"
    cdf = sc.zeros((N,2))
    cdf[:,0] = cdf_x; cdf[:,1] = cdf_y
    return cdf
    
def KS_stat(test_cdf, null_cdf):
    """ Gets the K-S statistic """
    D_max, D_sum, pos, j = -1.0, 0.0, -1, 0
    test_N, null_N = len(test_cdf), len(null_cdf)
    for i in range(test_N):
        # get difference if test data extends lower than null set
        if (test_cdf[i,0] < np.ma.min(null_cdf) ):
            D = test_cdf[i,1]  #test_y - 0.0, no value
        # get test difference if null data has non-zero value
        else:  #maybe a binary search would be better here?  Searching backwards from end
            for j in range(1, (null_N+1)):
                if null_cdf[(-1*j),0] <= test_cdf[i,0]:  break
            D = abs(test_cdf[i,1] - null_cdf[(-1*j),1])
        if D > D_max:  D_max = D; pos = [i, j]
        D_sum = D_sum + D
    print D_max, D_sum, pos
    return D_max, D_sum, pos

def kolmogorov_distribution():
    pass

def get_alpha(D_max, N):
    left = D_max*np.sqrt(N)
    # NEED TABLE OF VALUES!!!
    

def plot_cdfs(test_cdf, null_cdf, D_max, pos):
    # setup test cdf plot
    test_x, test_y = sc.zeros(2*len(test_cdf[:,0])+2), sc.zeros(2*len(test_cdf[:,0])+2)
    test_x[0], test_x[1] = (np.ma.min(test_cdf[:,0])-1.0), np.ma.min(test_cdf[:,0])
    test_y[0], test_y[1] = 0.0, 0.0
    for i in range(len(test_cdf[:,0])-1):
        test_x[(i*2)+2], test_x[(i*2)+3] = test_cdf[i,0], test_cdf[(i+1),0]
        test_y[(i*2)+2], test_y[(i*2)+3] = test_cdf[i,1], test_cdf[i,1]
    test_x[-2], test_y[-2] = test_cdf[-1,0], test_cdf[-1,1]
    test_x[-1], test_y[-1] = (test_cdf[-1,0]+1.0), test_cdf[-1,1]
    # setup null cdf plot
    null_x, null_y = sc.zeros(2*len(null_cdf[:,0])+2), sc.zeros(2*len(null_cdf[:,0])+2)
    null_x[0], null_x[1] = (np.ma.min(null_cdf[:,0])-1.0), np.ma.min(null_cdf[:,0])
    null_y[0], null_y[1] = 0.0, 0.0
    for i in range(len(null_cdf[:,0])-1):
        null_x[(i*2)+2], null_x[(i*2)+3] = null_cdf[i,0], null_cdf[(i+1),0]
        null_y[(i*2)+2], null_y[(i*2)+3] = null_cdf[i,1], null_cdf[i,1]
    null_x[-2], null_y[-2] = null_cdf[-1,0], null_cdf[-1,1]
    null_x[-1], null_y[-1] = (null_cdf[-1,0]+1.0), null_cdf[-1,1]
    # setup D plot
    D_x = [test_cdf[pos[0],0], test_cdf[pos[0],0]]
    D_y = [(test_cdf[pos[0],1]), (null_cdf[(-1*pos[1]),1])]
    # setup figure
    plt.figure(1)
    plt.plot(test_x, test_y, "k-")
    plt.plot(null_x, null_y, "k:")
    plt.plot(D_x, D_y, "r-")
    plt.xlabel(r"x")
    plt.ylabel("cumulative density")
    plt.show()
    plt.close("all")

if __name__ == "__main__":
    np.random.seed(int(ti.time()))
    dist_1 = get_CDF(np.random.normal(0.0, 20.0, 1000.0))
    dist_2 = get_CDF(np.random.normal(0.0, 20.0, 1000.0))
    max, sum, n = KS_stat(dist_1, dist_2)
    plot_cdfs(dist_1, dist_2, max, n)
    print "# --- Done"