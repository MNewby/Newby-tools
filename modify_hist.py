#! /usr/bin/env python  #'Bang' line - modify as needed

import math as m
import numpy as np
import scipy as sc
import hist as h
import files as f
import gctools_enc as gct
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt

"""Short script for modifying histgram data files.

Matthew Newby (RPI), March 22, 2011
"""

NAME = ['4147', '5024', '5053', '5272', '5466', '5904', '6205', '6341', '7078', '7089', 'Pal5']
DISTANCE = [19.3, 18.7, 18.5, 10.4, 15.6, 8.0, 7.7, 8.7, 11.0, 11.5, 21.0] 
AREA_IN = [0.0145, 0.1679, 0.0537, 0.1788, 0.0826, 0.076, 0.0664, 0.0175, 0.0504, 0.022, 0.0238]
AREA_BACK = [0.0915, 1.1548, 0.6922, 2.1823, 0.9737, 1.8537, 3.2279, 0.0635, 0.2517, 0.4120, 0.2740]
CUTOFF = [22.9, 22.9, 23.1, 22.5, 22.9, 21.7, 22.1, 0.0, 22.5, 22.7, 22.9]
CUT_INDEX = [39, 39, 40, 37, 39, 33, 35, 00, 37, 38, 39]
datafile = ['noU_NGC_4147_cluster.csv',
            'noU_NGC_5024_cluster.csv',
            'noU_NGC_5053_cluster.csv',
            'noU_NGC_5272_cluster.csv',
            'noU_NGC_5466_cluster.csv',
            'noU_NGC_5904_cluster.csv',
            'noU_NGC_6205_cluster.csv',
            'noU_NGC_6341_cluster.csv',
            'noU_NGC_7078_cluster.csv',
            'noU_NGC_7089_cluster.csv',
            'noU_Pal5_cluster.csv',
           ]
backfile = ['noU_NGC_4147_background.csv',
            'noU_NGC_5024_background.csv',
            'noU_NGC_5053_background.csv',
            'noU_NGC_5272_background.csv',
            'noU_NGC_5466_background.csv',
            'noU_NGC_5904_background.csv',
            'noU_NGC_6205_background.csv',
            'noU_NGC_6341_background.csv',
            'noU_NGC_7078_background.csv',
            'noU_NGC_7089_background.csv',
            'noU_Pal5_background.csv',
           ]
clus_wide = ['wide_4147_cluster.csv',
            'wide_5024_cluster.csv',
            'wide_5053_cluster.csv',
            'wide_5272_cluster.csv',
            'wide_5466_cluster.csv',
            'wide_5904_cluster.csv',
            'wide_6205_cluster.csv',
            'wide_6341_cluster.csv',
            'wide_7078_cluster.csv',
            'wide_7089_cluster.csv',
            'wide_Pal5_cluster.csv',
           ]
back_wide = ['wide_4147_background.csv',
            'wide_5024_background.csv',
            'wide_5053_background.csv',
            'wide_5272_background.csv',
            'wide_5466_background.csv',
            'wide_5904_background.csv',
            'wide_6205_background.csv',
            'wide_6341_background.csv',
            'wide_7078_background.csv',
            'wide_7089_background.csv',
            'wide_Pal5_background.csv',
           ]

from functions import SDSS_detection_efficiency
import gctools_enc as gce
### Get and plot ratio data
clusters = ['4147', '5024', '5053', 'Pal5']
styles = ['--', '-.', ':', '-']
data = []
for i in range(len(clusters)):
    data.append(f.read_data(clusters[i]+'_ratio_data_new.txt'))
plt.figure(1)
for i in range(len(clusters)):
    #plt.plot(data[i][:,0], data[i][:,1], label=clusters[i])
    plt.errorbar(data[i][:,0], data[i][:,1], yerr=data[i][:,2], label="NGC "+clusters[i], fmt=styles[i])
#Plot fit curve over data
x_func = sc.arange(19.0, 26.0, 0.1)
y_func = sc.zeros(len(x_func), float)
for i in range(len(x_func)):
    y_func[i] = gce.parabola_error(x_func[i])
x_sig = sc.arange(19.0, 26.0, 0.1)
y_sig = SDSS_detection_efficiency(x_sig)
plt.plot(x_func, y_func, 'k-', label='Functional Fit')
plt.plot(x_sig, y_sig, 'k:', label="Field Det. Eff.")
plt.xlabel(r'$g_0$', fontsize=14)
plt.ylabel('fractional completeness', fontsize=12)
plt.legend()
plt.show()

"""
### Compare first few steps of near avg and farther cluster; recaluculate farther cluster's
### bin heights to account for cutoff data at high magnitudes.
clus_name = 'Pal5'
avg_data = f.read_data('avg_normed_hist.txt')
far_data = f.read_data(clus_name+'_normed_hist.txt')
far_dist = DISTANCE[NAME.index(clus_name)]
modulus = (5.*(m.log10(far_dist*1000) - 1.))
x,y,z = [], [], []
for i in range(5):
    center = far_data[i,0] - modulus
    if round(center, 3) == round(avg_data[i,0], 3):
        x.append(center)
        y.append(far_data[i,1] - avg_data[i,1])
        z.append(far_data[i,1] / avg_data[i,1])
    else:  print "!!! No match!"
print x
print y
print z
print "# Average diff:", sc.mean(y)
print "# Diff std:", sc.std(y)
print "# Average Quot:", sc.mean(z)
print "# Quot std:", sc.std(z)
error_array = (np.sqrt( (avg_data[:,3])**2 + (far_data[:,2])**2 ) )
#print error_array
new_far = far_data[:,1] / sc.mean(z) #2.22447512508
#newer_far = far_data[:,1] - sc.mean(y)
#for i in range(len(newer_far)):
#    if newer_far[i] < 0.0:  newer_far[i] = 0.0
out_array = sc.zeros((len(far_data), 3), float)
out_array[:,0] = avg_data[:,0]
out_array[:,1] = new_far[:]
#h.plot_histogram(out_array, 0.2, limits = [], name='far_clus_renormed', x_label=r'$M_g$')
h.plot_multiple_hist(out_array[:,0], [avg_data[:,1], new_far], 0.2,
                     limits = [], name=(clus_name+'_and_avg_two'), x_label=r'$M_g$')
ratio_data = []
for i in range(len(new_far)):
    ratio_data.append( [ far_data[i,0], (new_far[i] / avg_data[i,1]), error_array[i] ])
f.write_data(sc.array(ratio_data), (clus_name+"_ratio_data_new.txt") )
"""

"""
### Get histogram files in, convert to absolute magnitude, then average
file_names = ['6205_normed_hist.txt', '5904_normed_hist.txt', '5272_normed_hist.txt']
dist = [7.7, 8.0, 10.4]
counts = [2164, 3698, 4364]
data = []
for cluster in file_names:
    data.append(f.read_data(cluster))
# Get absolute magnitudes and check that they are equal
abs_mags = []
for i in range(len(data)):
    modulus = (5.*(m.log10(dist[i]*1000) - 1.))
    abs_mags.append(data[i][:,0] - modulus)
for i in range(len(abs_mags[0])):
    if (round(abs_mags[0][i], 3) == round(abs_mags[1][i], 3) == round(abs_mags[2][i], 3)):
        print "#--Absolute magnitudes check out"
    else:
        print "!!! Absolute magnitudes not good!"
        print abs_mags[0][i], abs_mags[1][i], abs_mags[2][i]
# Get counts and average them
avg_list = []
for i in range(len(data[0][:,1])):
    set = [data[0][i,1], data[1][i,1], data[2][i,1]]
    avg_counts = sc.mean(set)
    avg_std = sc.std(set)
    err_set = [(set[0]*counts[0]), (set[1]*counts[1]), (set[2]*counts[2])]
    for j in range(len(err_set)):
        err_set[j] = (np.sqrt(err_set[j] + 1.0) + 1.0) / counts[j]
    print err_set
    avg_error = np.sqrt( (err_set[0])**2 + (err_set[1])**2 + (err_set[2])**2 ) / 3.0
    print avg_error
    avg_list.append([avg_counts, avg_std, avg_error])
avg_array = sc.array(avg_list)
out_array = sc.zeros((len(abs_mags[0]),4), float)
out_array[:,0] = abs_mags[0]
out_array[:,1] = avg_array[:,0]
out_array[:,2] = avg_array[:,1]
out_array[:,3] = avg_array[:,2]
f.write_data(out_array, fileout=('avg_normed_hist.txt'), delimiter='\t',
             header=("bin center (abs), average normed counts, std_counts, sigma"))
h.plot_histogram(out_array, 0.2, limits = [], name='avg_clus_normed', x_label=r'$M_g$')
"""

"""
### Get distribution of stars within two magnitudes for one cluster
cluster = 'Pal5'
# Get data and absolute magnitude limits for a cluster, given it's name and distance
index = NAME.index(cluster)
clus_data = f.read_csv(clus_wide[index])
abs_limit = [3.5675463741375903, 7.5675463741375903]
modulus = (5.*(m.log10(DISTANCE[index]*1000) - 1.))
app_limit = [(abs_limit[0] + modulus), (abs_limit[1] + modulus)]
print "#---apparent magnitude limits:", app_limit
clus_list = []
for i in range(len(clus_data[:,0])):
    g_mag = clus_data[i,2]
    g_minus_r = (clus_data[i,2] - clus_data[i,3])
    u_minus_g = (clus_data[i,4] - clus_data[i,2])
    if u_minus_g < 0.4:  continue
    if g_minus_r < 0.0:  continue
    if g_minus_r > 0.8:  continue
    if g_mag > app_limit[1]: continue
    if g_mag < app_limit[0]: continue
    abs_mag = g_mag - (5.*(m.log10(DISTANCE[index]*1000) - 1.))
    clus_list.append([g_mag, abs_mag, g_minus_r, u_minus_g])
total = len(clus_list)
clus_array = sc.array(clus_list)
clus_hist = h.make_hist(clus_array[:,0], 0.2, app_limit)
#for i in range(len(clus_hist[:,0])):
#    clus_hist[i,1] = (clus_hist[i,1] / total)
data_out = sc.zeros((len(clus_hist[:,0]), 3), float)
for i in range(len(clus_hist[:,0])):
    data_out[i,0] = clus_hist[i,0]
    data_out[i,1] = (clus_hist[i,1] / total)
    data_out[i,2] = ((np.sqrt(clus_hist[i,1]+1.0) + 1.0) / total)
f.write_data(data_out, fileout=(cluster+'_normed_hist.txt'), delimiter='\t', header=(cluster+' stars: '+str(total)))
h.plot_histogram(clus_hist, 0.2, app_limit, name=cluster+'_clus_normed_hist', x_label='g')
"""


'''### Get distributions of stars between background and cluster, for stars
### between 22.0 < g < 25.0 and 0.1 < g-r < 0.3
# Only want:
##5272 - 3
#5904 - 5
#6205 - 6
###7078 - 8
#7089 - 9
bin_size = 0.2
runs = []
for i in range(len(NAME)):
    if i < 3:  continue 
    if i == 4:  continue
    if i == 7:  continue
    if i > 9: continue
    factor = AREA_IN[i] / AREA_BACK[i]
    clus_data = f.read_csv(datafile[i])
    back_data = f.read_csv(backfile[i])
    clus_list = []
    for j in range(len(clus_data[:,0])):
        if (clus_data[j,2] < 22.0) or (clus_data[j,2] > 25.0):  continue
        g_minus_r = (clus_data[j,2] - clus_data[j,3])
        if (g_minus_r < 0.1) or (g_minus_r > 0.3):  continue
        clus_list.append([clus_data[j,2], g_minus_r])
    back_list = []
    for j in range(len(back_data[:,0])):
        if (back_data[j,2] < 22.0) or (back_data[j,2] > 25.0):  continue
        g_minus_r = (back_data[j,2] - back_data[j,3])
        if (g_minus_r < 0.1) or (g_minus_r > 0.3):  continue
        back_list.append([back_data[j,2], g_minus_r])
    clus_hist = h.make_hist(sc.array(clus_list)[:,0], bin_size, spread=[22.0, 25.0])
    back_hist = h.make_hist(sc.array(back_list)[:,0], bin_size, spread=[22.0, 25.0])
    new_hist = sc.zeros((clus_hist.shape), float)
    for j in range(len(clus_hist[:,0])):
        #print clus_hist[j,1], back_hist[j,1]
        if (clus_hist[j,0] == back_hist[j,0]):
            new_hist[j,0] = clus_hist[j,0]
            if (back_hist[j,1] > 0.0):  
                #new_hist[j,1] = (clus_hist[j,1] / (factor*back_hist[j,1]))
                new_hist[j,1] = (clus_hist[j,1] / back_hist[j,1])
            else:
                new_hist[j,1] = 0.0
        else: print 'OH NO!'
    #print new_hist
    runs.append(sc.array(new_hist))
    #h.plot_histogram(clus_hist, bin_size, limits=[22.0, 25.0], name=(NAME[i]+'_clus_hist_out'), x_label='g')
    #h.plot_histogram(back_hist, bin_size, limits=[22.0, 25.0], name=(NAME[i]+'_back_hist_out'), x_label='g')
    #h.plot_histogram(new_hist, bin_size, limits=[22.0, 25.0], name=(NAME[i]+'_new_nofac'), x_label='g')
plt.figure(1)
for i in range(len(runs)):
    plt.plot(runs[i][:,0], runs[i][:,1])
plt.xlim(22.0, 25.0)
plt.show()
'''

''' ### Plot representations of histograms, of all clusters in same plot
magnitudes, percentages, counter = [], [], 0
for i in range(len(NAME)):
    if i == 7: continue  #skip 6341
    clus_file = NAME[i]+'_back_cumulative.txt'
    clus_data = f.read_data(clus_file)
    magnitudes.append(clus_data[:,0])
    percentages.append(clus_data[:,1])
    #print percentages
    #total = sc.sum(percentages[counter])
    total = percentages[counter][-1]
    print total
    for j in range(len(percentages[counter])):
        percentages[counter][j] = (percentages[counter][j]/total)
    counter = counter + 1
counter = 0
plt.figure(1)
for i in range(len(NAME)):
    if i == 7: continue  #skip 6341
    plt.plot(magnitudes[counter], percentages[counter])
    #x = DISTANCE[i]*sc.ones(len(magnitudes[counter]))
    #y = magnitudes[counter]
    #z = 1000*percentages[counter]
    #plt.scatter(x[:CUT_INDEX[i]], y[:CUT_INDEX[i]], z[:CUT_INDEX[i]], 'k', 'o')
    #plt.scatter(x, y, z, 'k', 'o')
    counter = counter + 1
plt.show()
'''
    

'''  ###  Get magnitude that contains 95% of the population
for i in range(len(NAME)):
    if i == 7: continue  #skip 6341
    clus_file = NAME[i]+'_back_cumulative.txt'
    clus_data = f.read_data(clus_file)
    total = clus_data[-1,1]
    for j in range(len(clus_data[:,0])):
        percent = (clus_data[j,1] / total)
        if percent > 0.95:
            print NAME[i], percent, clus_data[j,0], DISTANCE[i], j, len(clus_data[:,0])
            break
'''

''' ### Subtract background from cluster
for i in range(len(NAME)):
    if i == 7: continue  #skip 6341
    factor = AREA_IN[i] / AREA_BACK[i]
    clus_file = NAME[i]+'_normal.txt'
    back_file = NAME[i]+'_back_normal.txt'
    clus_data = f.read_data(clus_file)
    back_data = f.read_data(back_file)
    # New histograming program!
    new_hist = sc.zeros(clus_data.shape)
    new_hist[:,0] = clus_data[:,0]
    for j in range(len(clus_data[:,0])):
        new_hist[j,1] = clus_data[j,1] - (factor*back_data[j,1])
        if (new_hist[j,1] < 0.0):  new_hist[j,1] = 0.0
    h.plot_histogram(new_hist, 0.2, limits=[], name=(NAME[i]+'_combined'), x_label=r'$g$')
    f.write_data(new_hist, fileout=(NAME[i]+'_combined.txt') )
    print '#-Modulus for', NAME[i], (5.*(m.log10(DISTANCE[i]*1000) - 1.))
'''

"""  ###  Testing different histogram methods for consistancy; appears correct. -March 23, 2011
datafile = 'noU_NGC_6205_cluster.csv'
backfile = 'noU_NGC_6205_background.csv'
clus_data = f.read_csv(datafile)
back_data = f.read_csv(backfile)
x_clus = clus_data[:,2]
x_back = back_data[:,2]

bin_size = 0.2
i = 6  #NGC 6205
factor = AREA_IN[i] / AREA_BACK[i]

modulus = (5.*(m.log10(DISTANCE[i]*1000) - 1.))
clus_hist = h.make_hist(x_clus, bin_size, spread=[15.0, 25.0])
back_hist = h.make_hist(x_back, bin_size, spread=[15.0, 25.0])
mod_hist = sc.zeros(clus_hist.shape)
for j in range(len(clus_hist[:,0])):
    if clus_hist[j,0] != back_hist[j,0]:  print  '!!! NOT EQUAL!!!', clus_hist[j,0], back_hist[j,0]
    mod_hist[j,0] = clus_hist[j,0]
    mod_hist[j,1] = clus_hist[j,1] - (factor*clus_hist[j,1])
    if (mod_hist[j,1] < 0.0):  mod_hist[j,1] = 0.0
h.plot_histogram(clus_hist, bin_size, limits=[15.0, 25.0], name=(NAME[i]+'_clus'), x_label=r'$g$')
h.plot_histogram(back_hist, bin_size, limits=[15.0, 25.0], name=(NAME[i]+'_back'), x_label=r'$g$')
h.plot_histogram(mod_hist, bin_size, limits=[15.0, 25.0], name=(NAME[i]+'_combined'), x_label=r'$g$')

old_hist = gct.binned_background(x_clus, bin_size=bin_size, cut_hist=0, in_area=AREA_IN[i],
                                 back_area=AREA_BACK[i], low_cut=15.0, high_cut=25.0)
#oback_hist = gct.binned_background(x_back, bin_size=bin_size, cut_hist=0, in_area=AREA_IN[i],
#                                 back_area=AREA_BACK[i], low_cut=15.0, high_cut=25.0)
mix_hist = gct.binned_background(x_clus, x_back, bin_size, cut_hist=0, in_area=AREA_IN[i],
                                 back_area=AREA_BACK[i], low_cut=15.0, high_cut=25.0)
h.plot_histogram(old_hist, bin_size, limits=[15.0, 25.0], name=(NAME[i]+'_oriclus'), x_label=r'$g$')
h.plot_histogram(mix_hist, bin_size, limits=[15.0, 25.0], name=(NAME[i]+'_oricom'), x_label=r'$g$')
"""

print '#---Done'