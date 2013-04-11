import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func

'''python script for acquiring only 2-sigma significant points in a CMD.
Matthew Newby, August 2, 2011'''


def rejection(filename, min=None, max=None, bin=0.2, reject=2.0, thresh=1000.0):
    data = fi.read_data(filename, ",")
    if min == None:
        min = np.ma.min(data[:,2])
    if max == None:
        max = np.ma.max(data[:,2])
    data_out = []
    steps = int((max-min)/bin)+1
    for i in range(steps):
        # build a data set for each bin
        data_temp = []
        min_temp, max_temp = (min + i*bin), (min + (i+1)*bin)
        for j in range(len(data[:,0])):
            if (data[j,2] >= min_temp) and (data[j,2] < max_temp):
                data_temp.append(data[j,:])
        print "number of stars in bin {0}: {1}".format(i, len(data_temp))
        if len(data_temp)==0:  continue  #skip empty sets
        # continue checking for rejections until no rejections occur
        kills = [1] 
        while len(kills) > 0:
            kills, new_temp = [], []
            # pack the g-r values into a single list
            g_minus_r = []
            for j in range(len(data_temp)):
                g_minus_r.append(data_temp[j][2]-data_temp[j][3])
            g_minus_r = sc.array(g_minus_r)
            avg = sc.mean(g_minus_r)
            stdev = sc.std(g_minus_r)
            #print avg, stdev
            # Find the values that lie outside 'reject' sigma
            for j in range(len(g_minus_r)):
                diff = ma.fabs( g_minus_r[j] - avg )
                if diff > (reject*stdev):  kills.append(data_temp[j])
                elif diff > thresh:  kills.append(data_temp[j])
                else:  new_temp.append(data_temp[j])
            data_temp = new_temp
        data_out = data_out + data_temp
        print "Stars kept in bin: {0}".format(len(data_temp))
    # restructure data
    data_out = sc.array(data_out)
    print " {0} stars remaining out of {1} original stars".format(len(data_out[:,0]), len(data[:,0]))
    g_minus_r_ori = data[:,2] - data[:,3]
    g_minus_r_out = data_out[:,2] - data_out[:,3]
    # show plot, then ask to save data set
    fig = plt.figure()
    #plt.scatter(g_minus_r_ori, data[:,2], 2, marker='o', facecolor="r", edgecolor='face')
    plt.scatter(g_minus_r_out, data_out[:,2], 1, 'k', 'o')
    plt.ylim(max, min)
    plt.xlim(-1.0, 2.0)
    plt.ylabel(r'$g_0$')
    plt.xlabel(r"$(g-r)_0$")
    plt.show()
    fname = raw_input("save data? [filename to save, empty to skip]")
    if fname != "":
        fi.write_data(data_out, fname, header="#")
    print "# --- Done"
    
if __name__ == "__main__":
    rejection("NGC_6205_all_cluster.csv", 15.0, 26.0, bin=0.2, reject=2.5, thresh=0.5)