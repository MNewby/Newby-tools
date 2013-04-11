#! /usr/bin/env python

import files as f
import scipy as sc
import math as m
import numpy as np
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt


"""This program plots data sets from isocrones in a pretty ganged plot.

Matthew Newby (RPI), Oct 2, 2010
"""

#Could roll things up in lists.
save = 1
save_name = 'Iso_combined_fit_noU_final.ps'
gc_names = ['4147', '5024', '5053', '5272', '5466', '5904', '6205', '6341', '7078', '7089', 'Pal 5']
gc_distances = [19.3, 18.7, 18.5, 10.4, 15.6, 8.0, 7.7, 8.7, 11.0, 11.5, 21.0]
#gc_metals = [-1.58, -1.89, -2.71, -1.43, -2.14, -1.17, -1.42, -2.17, -2.04, -1.38, -1.24]
fid_res = [5, 2, 2, 2, 2, 2, 2, 5, 2, 2, 2]
cutoff = [3.0, 3.2, 3.1, 3.2, 3.2, 3.2, 3.3, 3.2, 3.0, 3.3, 3.4]

pos_x = 1
pos_y = 1
subs = []
"""Plot Initializations"""
fig = plt.figure()
plt.subplots_adjust(hspace=0.001, wspace=0.001)
#plt.title(save_name)
#plt.xlabel('g-r')
#plt.ylabel('Mg')
for i in range(len(gc_names)):
    name = gc_names[i]
    distance = gc_distances[i]
    #metal = gc_metals[i]
    HB_limit = cutoff[i]
    
    """Get Data"""
    if (i == 10):
        data_iso = f.read_data( ('Iso_new_Pal5.dat') )
        data_clus = f.read_csv('noU_Pal5_cluster.csv')  #'HR_Pal_5_cluster.csv'
        data_fid = f.read_data( ('Pal5_2_fiducial_out.txt') )
    else:
        data_iso = f.read_data( ('Iso_new_'+name+'.dat') )
        #data_iso = f.read_data( ('Iso_series_'+name+'_A.dat') )
        data_clus = f.read_csv('noU_NGC_'+name+'_cluster.csv') #'HR_NGC_'+name+'_cluster.csv'
        data_fid = f.read_data( ('NGC_'+name+'_'+str(fid_res[i])+'_fiducial_out.txt') )

    """Setup Data for Plotting"""
    mm = -0.01463023
    bb = 0.08928602
    iso_x_list, iso_y_list = [], []
    for j in range(len(data_iso[:,8])):
        if (data_iso[j,8] > HB_limit):
            iso_x_list.append((data_iso[j,8] - data_iso[j,9]) + mm*data_iso[j,8] + bb)
            iso_y_list.append( data_iso[j,8] )
    iso_x = sc.array(iso_x_list)
    iso_y = sc.array(iso_y_list)
    clus_x = data_clus[:,2] - data_clus[:,3]
    clus_y = data_clus[:,2] - 5.*(m.log10(distance*1000) - 1.0) #'g' values
    fid_x_list, fid_y_list = [], []
    for j in range(len(data_fid[:,0])):
        if (data_fid[j,0] > HB_limit):
            fid_x_list.append(data_fid[j,1])
            fid_y_list.append(data_fid[j,0])
    fid_x = sc.array(fid_x_list)
    fid_y = sc.array(fid_y_list)

    """Actual Plotting Code"""
    yplots = [0,4,8]  #plots that will contain the y axes, -1 to give correct index
    #print pos_x, pos_y, i
    if (i > 0):
        if (i == 4 or i == 8):  sp = plt.subplot(3,4,(i+1), sharex=subs[(pos_x-1)])
        elif (pos_y > 1):  sp = plt.subplot(3,4,(i+1), sharex=subs[(pos_x-1)] ,
                                          sharey=subs[yplots[(pos_y-1)]])
        else:  sp = plt.subplot(3,4,(i+1), sharey=subs[0])
    else:
        sp = plt.subplot(3,4,(i+1))
    subs.append(sp)
    #plt.errorbar(data_x, data_y, yerr=y_err, xerr=x_err, ecolor='b', fmt=None)
    #plt.plot(np.sort(data_x), np.sort(line_y), c='g')
    sp.scatter(clus_x, clus_y, c='green', s=0.1, marker='o')
    sp.plot(iso_x, iso_y, 'b-')
    sp.scatter(fid_x, fid_y, c='red')
    plt.xticks(sc.arange(0.1, 0.6, 0.1), fontsize=6)
    plt.yticks(sc.arange(2.0, 8.0, 1.0), fontsize=6)
    plt.ylim(8.0, 1.0)  #Maybe need to change these to 'sp'?
    plt.xlim(0.0, 0.6)
    """ Setup labels """
    if pos_y < 3:
        if (i!=7): plt.setp(sp.get_xticklabels(), visible=False)
    if (pos_y == 3 or i==7):  plt.xlabel('$g_0-r_0$', fontsize=6)
    if pos_x > 1:  plt.setp(sp.get_yticklabels(), visible=False)
    if (pos_x == 1):  plt.ylabel('$M_g$', fontsize=6)
    if (i == 10):  plt.text(0.1, 2.0, name, fontsize=6)
    else:  plt.text(0.05, 2.0, ('NGC '+name), fontsize=6)

    """Keeps Track of which Axes to Copy """
    pos_x = pos_x + 1
    if pos_x > 4:
        pos_x = 1
        pos_y = pos_y + 1
    if pos_y > 3:  print '!!!---Error! pos_y greater than 3!'

"""Output"""
if save == 1:
    plt.savefig(save_name, papertype='letter')
    print '#---data plot saved as', save_name
else:
    plt.show()
plt.close('all')
print '#---Done with data plot'