#! /usr/bin/env python

import csv
import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt

"""This program contains code for creating fiducial sequences from color-
magnitude-diagrams.

Based on the 'Isocrone_fit.py' script
Matthew Newby (RPI), Oct 11, 2010
"""

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
data_file = 'HR_Pal_5_cluster.csv'
cluster_name = 'Pal5_2'
distance = 21.0  #starting distance
step_size = 0.2  #smaller gives more points, but larger errors
save_data = 1  #saves isocrone files as output if set to 1
#initializations
plot_title = cluster_name + ', fiducial fit(green)'
file_str = 'HR_fiducial_' + cluster_name + '.ps' #plot file string
#load in and prepare SDSS data
gc_data = f.read_csv(data_file)
gc_l, gc_w = gc_data.shape
gc_array = sc.zeros((gc_l,2), float)
gc_array[:,0] = gc_data[:,2] - 5.*(m.log10(distance*1000) - 1.0) #'g' values  
gc_array[:,1] = (gc_data[:,2] - gc_data[:,3])  #'g-r' values
bounds = sc.arange(1.0, 7.0, step_size)
#break up data CMD into strips in Mg
#reduce Mg strips - reject outliers, get mean, st_dev for each strip
#Note - if distances are fitted, then this part needs to be done every iteration.
fit_array = sc.zeros((len(bounds), 3),float)  #average g, g-r mean, g-r st_dev
for i in range(1, len(bounds)):
    star_list, g_list = [], []
    for j in range(gc_l):
        if ( (gc_array[j,0] < bounds[i]) and (gc_array[j,0] >= bounds[i-1])):
            #Do standard deviation in g, too?
            g_list.append(gc_array[j,0])
            star_list.append(gc_array[j,1])
    #3 sigma rejection
    reject = 1
    while (reject == 1):
        if (star_list == []):
            print 'empty star list!!!  between:', bounds[i], bounds[i-1], i
            fit_array[(i-1),:] = fit_array[(i-2),:]
            break
        if (len(star_list) == 1):
            print 'single star in star list!!!'
            fit_array[(i-1),:] = g_list[0], star_list[0], 0.1
            break
        fit_array[(i-1),0] = get_mean(g_list)
        fit_array[(i-1),1] = get_mean(star_list)
        fit_array[(i-1),2] = get_st_dev(star_list, fit_array[(i-1),1])  
        kill_list = []
        for j in range(len(star_list)):
        #now change errors from distribution widths to uncertainty in means
            if (abs(star_list[j] - fit_array[(i-1),1]) > (2.0*fit_array[(i-1),2]) ):
                kill_list.append(j)
        if (kill_list == []):
            fit_array[(i-1),2] = (fit_array[(i-1),2] / float(len(star_list)))
            reject = 0
        else:
            print 'kill length', len(kill_list), 'bin', (i)
            for j in range(len(kill_list)):
                star_list.pop(kill_list[j])
                g_list.pop(kill_list[j])
                for k in range(len(kill_list)):  kill_list[k] = (kill_list[k]-1)
#now have 'iso_array' and 'fit_array' - save to file!
if (save_data==1):
    outname = cluster_name + '_fiducial_out.txt'
    dist_str = (str(distance))
    header_txt = 'g, g-r, uncertainty mean, from file ' + data_file + ' with distance ' + dist_str
    #print out_data
    #print out_data.shape
    if (f.write_data(fit_array, outname, header=header_txt) == 1):
        print 'iso g, g-r, sigma, data successfully saved'
#fit distances?  Hard!  Need to start from beginning..."""
plt.figure()
y = gc_array[:,0]
x = gc_array[:,1]
plt.scatter(x,y, c='black', marker='+', label='HR data')
#plt.errorbar(x,y,xerr=fit_array[:,2])
plt.errorbar(fit_array[:,1],fit_array[:,0],xerr=fit_array[:,2], fmt=None, ecolor='g')
plt.scatter(fit_array[:,1],fit_array[:,0], c='green', label='fiducial sequence')
plt.title(plot_title, fontsize='small')
plt.xlabel('g-r')
plt.ylabel('Mg')
leg = plt.legend(loc='best')
for t in leg.get_texts():
    t.set_fontsize('small')
plt.ylim(7.0,1.0)
plt.xlim(0.0, 0.6)
plt.savefig(file_str, papertype='letter')
#plt.show()
plt.close('all')
print '-Done'