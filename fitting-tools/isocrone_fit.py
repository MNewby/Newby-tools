#! /usr/bin/env python

import csv
import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt

"""This program contains code for building a grey atmospere 'fudge function'
into the Girardi isocrones for ugriz filters.  Use this one.

Matthew Newby (RPI), Aug 20, 2010
"""
#plot './stars_ngc_6205.txt' using ($3-$4):3, './6205_isocrone_12.3gyr_1.41feH.dat' using ($9-$10):9

def get_mean(in_list):
    mean = 0.0
    for j in range(len(in_list)):
        mean = mean + in_list[j]
    return ( mean / float((len(star_list))) )
    
def get_st_dev(in_list, mu):
    st_dev = 0.0
    for j in range(len(in_list)):
        st_dev = st_dev + ((in_list[j] - mu)**2)
    return np.sqrt(st_dev/float(len(in_list)))
    
#initial parameters; this would be a good point for a wrapper.
iso_file = '6205_isocrone_12.3gyr_1.41feH.dat'
data_file = 'HR_NGC_6205_cluster.csv'
cluster_name = 'NGC_6205_CG97_clean2'
distance = 7.7  #starting distance
save_data = 0  #saves isocrone files as output if set to 1
#initializations
#plot_title = cluster_name + ', isocrone(blue) and fiducial fit(green)'
file_str = 'HRISO_compare_' + cluster_name  #plot file string
#load in and prepare SDSS data
gc_data = f.read_csv(data_file)
gc_l, gc_w = gc_data.shape
gc_array = sc.zeros((gc_l,2), float)
gc_array[:,0] = gc_data[:,2] - 5.*(m.log10(distance*1000) - 1.0) #'g' values  
gc_array[:,1] = (gc_data[:,2] - gc_data[:,3])  #'g-r' values
#load in and prepare isocrone data
iso_data = f.read_data(iso_file)  #<--might fail due to leading whitespace?
iso_l, iso_w = iso_data.shape
iso_in = sc.zeros((iso_l,2), float)
iso_in[:,0] = iso_data[:,8]  #'Mg' values
iso_in[:,1] = (iso_data[:,8] - iso_data[:,9])  #'g-r' values
#Chop up isocrone so that the bounds are equal
gminusr = 0.6
for i in range(1,iso_l):
    if ( (iso_in[i,1] < gminusr) and (iso_in[(i-1),1] > gminusr) ):
        high_g = i
    if ( (iso_in[i,1] > gminusr) and (iso_in[(i-1),1] < gminusr) ):
        low_g = i
        break
iso_array = iso_in[high_g:low_g,:]
iso_l, iso_w = iso_array.shape  
#break up data CMD into strips in Mg
bounds = []
for i in range(iso_l):
    if (i == 0):
        bounds.append((iso_array[0,0] + abs(iso_array[0,0]-iso_array[0,1])))
        continue
    bounds.append(iso_array[i,0] + (abs(iso_array[i-1,0] - iso_array[i,0])/2.0))
    if (i == (iso_l - 1)):
        bounds.append((iso_array[-1,0] - abs(iso_array[-1,0]-iso_array[-2,0])))
        continue
#reduce Mg strips - reject outliers, get mean, st_dev for each strip
#Note - if distances are fitted, then this part needs to be done every iteration.
fit_array = sc.zeros((iso_l, 3),float)  #average g, g-r mean, g-r st_dev
for i in range(1, len(bounds)):
    star_list, g_list = [], []
    for j in range(gc_l):
        if ( (gc_array[j,0] < bounds[i-1]) and (gc_array[j,0] >= bounds[i])):
            #Do standard deviation in g, too?
            g_list.append(gc_array[j,0])
            star_list.append(gc_array[j,1])
    #3 sigma rejection
    reject = 1
    while (reject == 1):
        if (star_list == []):
            print 'empty star list!!!  between:', bounds[i-1], bounds[i], i
            fit_array[(i-1),:] = fit_array[(i-2),:]
            break
        if (len(star_list) == 1):
            print 'single star in star list!!!'
            fit_array[(i-1),:] = g_list[0], star_list[0], 0.1
            break
        fit_array[(i-1),0] = get_mean(g_list)
        fit_array[(i-1),1] = get_mean(star_list)
        fit_array[(i-1),2] = get_st_dev(star_list, fit_array[(i-1),1])  #good to here----
        kill_list = []
        for j in range(len(star_list)):
        #now change errors from distribution widths to uncertainty in means
            if (abs(star_list[j] - fit_array[(i-1),1]) > (2.0*fit_array[(i-1),2]) ):
                kill_list.append(j)
        if (kill_list == []):
            fit_array[(i-1),2] = (fit_array[(i-1),2] / float(len(star_list)))
            reject = 0
        else:
            print 'kill length', len(kill_list), 'bin', (i-1)
            for j in range(len(kill_list)):
                star_list.pop(kill_list[j])
                g_list.pop(kill_list[j])
                for k in range(len(kill_list)):  kill_list[k] = (kill_list[k]-1)
#Fit strips to Girardi isocrones using different fit functions
#now have 'iso_array' and 'fit_array' - save to file!
if (save_data==1):
    outname = cluster_name + '_out.txt'
    out_data = sc.zeros((iso_l, 4),float)
    out_data[:,0], out_data[:,1], out_data[:,2], out_data[:,3] = \
        iso_array[:,0], iso_array[:,1], fit_array[:,1], fit_array[:,2]
    dist_str = (str(distance))
    header_txt = 'g, iso g-r, fiducial g-r, uncertainty mean, from file ' + iso_file + ' with distance ' + dist_str
    #print out_data
    #print out_data.shape
    if (f.write_data(out_data, outname, header=header_txt) == 1):
        print 'iso g, iso g-r, fid g-r, fiducial sigma, data successfully saved'
#fit distances?  Hard!  Need to start from beginning..."""
plt.figure()
A = -0.01463
B = 0.08929
y = iso_array[:,0]
x = (fit_array[:,1] - iso_array[:,1])
plt.scatter(x,y, c='black', marker='+', label='Subtracted data')
plt.plot([0.0,0.0], [np.ma.min(iso_array[:,0]),np.ma.max(iso_array[:,0])], 'r-', label='_nolegend_')
plt.plot([(A*np.ma.min(y)+B),(A*np.ma.max(y)+B)], [np.ma.min(y), np.ma.max(y)], 'g:', label='Line Fit')
#plt.plot([(A*np.ma.min(y)+B),np.ma.min(y)], [(A*np.ma.max(y)+B), np.ma.max(y)], 'g:', label='Line Fit')
#plt.errorbar(x,y,xerr=fit_array[:,2])
plt.errorbar(iso_array[:,1],iso_array[:,0], fmt='o', mfc='c', mec='c', mew=1, label='Girardi isocrone')
plt.errorbar(fit_array[:,1],iso_array[:,0],xerr=fit_array[:,2], fmt=None, ecolor='k')
plt.scatter(fit_array[:,1],iso_array[:,0], marker='d', c='black', label='Fiducial sequence')
#plt.title(plot_title, fontsize='small')
plt.ylim(np.ma.max(iso_array[:,0]), np.ma.min(iso_array[:,0]))
plt.xlabel('$g_0 - r_0$')
plt.ylabel('$M_g$')
#locs, lbls = plt.xticks()
#plt.xticks(locs, fontsize=8)
#locs, lbls = plt.yticks()
#plt.yticks(locs, fontsize=8)
leg = plt.legend(loc='lower center')
for t in leg.get_texts():
    t.set_fontsize(8)
#plt.savefig(file_str, papertype='letter')
plt.show()
plt.close('all')
print '-Done'