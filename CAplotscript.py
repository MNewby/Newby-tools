import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt

'''python script for plotting convolved analysis data
Matthew Newby, August 9, 2010'''

clusters = ['4147', '5024', '5053', '5272', '5466', '5904', '6205', '7078', '7089', 'Pal5']
file_str = '_CA_Harris.txt'
plot_str = '_convolve_analyzed.PS'
for i in range(len(clusters)):
    print '#-working on cluster', clusters[i]
    filename = clusters[i] + file_str
    raw_data = f.read_data(filename)
    l, w = raw_data.shape
    con_data, flux_data = sc.zeros( ((l/3), w)), sc.zeros( ((l/3), w))
    #print con_data.shape, flux_data.shape
    x, y, z = 1, 0, 0
    for j in range(l):
        if (x==3):
            #print j, x, l, z
            flux_data[z,:] = raw_data[j,:]
            x, z = 1, (z+1)
            continue
        if (x==2):
            #print j, x, y
            con_data[y,:] = raw_data[j,:]
            x, y = (x+1), (y+1)
            continue
        if (x==1):  x = (x + 1)
    #print con_data
    #print flux_data
    plt.figure(1)
    title_str = 'plot of ' + clusters[i]
    plt.title(title_str)
    plt.subplot(221)
    plt.title('g-r < 0.0')
    plt.scatter(con_data[:,0], con_data[:,1])
    plt.subplot(222)
    plt.title('f-turnoff')
    plt.scatter(con_data[:,0], con_data[:,2])
    plt.subplot(223)
    plt.title('g-r > 0.6')
    plt.scatter(con_data[:,0], con_data[:,3])
    plt.subplot(224)
    plt.title('Flux, from left=square. f-turn-off=circle, from right=triangle')
    plt.scatter(flux_data[:,0], flux_data[:,1], c='b', marker='s')
    plt.scatter(flux_data[:,0], flux_data[:,2], c='g', marker='o')
    plt.scatter(flux_data[:,0], flux_data[:,3], c='r', marker='^')
    save_str = clusters[i] + plot_str
    plt.savefig(save_str, papertype='letter')
    plt.close('all')
    print '#-plot successfully saved'
print '#-all jobs complete'