import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import sys

'''python script for plotting gc data; this is the one to use!!!
Matthew Newby, August 1, 2010

Takes 2 argv:  in file and an identifier
'''

def sigmoid_function(x, params):
        #Looks good, takes three parameters
        data = ( params[0]/(1 + np.exp(-1.0*(x-params[1])) ) ) + params[2]
        return data
def linear_function(x, params):
        #looks like this function is ok
        y = (x*params[0] + params[1])
        return y
def const_linear_function(x, params):
        #looks like this function is ok
        #y = ((x-20.0)*params[0] + 1.24127)
        y = ((x-18.0)*params[0] + 1.24009)
        return y


""" Initial Inputs """
results_file = sys.argv[1]
identifier = sys.argv[2]
xoff = -1.0  #text position offsets
yoff = 0.01
find_average = 0  #change
file_to_average = results_file  #This is the file to be used in finding averages
plot_errors = 0  #use 1 to plot error bars
# Loop through mu, sig_l, sig_r, plots
two_means = [0, 0, 4.31356, 0.4199, 0]
for plot in range(2,5):
    x = 0       #this uses distance by default
    y = plot	#change always
    e = y + 4
    #File names for input and output
    #results_file = 'mid_bin_results.txt'
    if (y==2):  leader = 'mu_'
    elif (y==3):  leader = 'sigl_'
    elif (y==4):  leader = 'sigr_'
    else:  leader = 'something_'
    plot_file = leader+identifier+'_GCs.ps'
    if (find_average == 1):
        averages_file = leader+identifier+'_GCs.txt'
    #Initializing variables
    NAME = ['4147', '5024', '5053', '5272', '5466', '5904', '6205', '7078', '7089', 'Pal5'] #'6341',
    #NAME = ['5272', '5466', '5904', '6205', '7078', '7089'] # Short Series
    axes_names = ['$d_{eff}$ (kpc)', 'Fwhm', '$\mu$', '$\sigma_l$', '$\sigma_r$',
                  'Amplitude', '$\mu$ Error', '$\sigma_l$ Error', '$\sigma_r$ Error', 'Amplitude Error']
    formats = ['b:', 'g-', 'r:', 'k-', 'c-', 'm-', 'b-', 'r:', 'g-', 'k-']
    #formats = ['b-', 'g-', 'r-', 'c-', 'r:', 'p-'] # Short series
    x_label = axes_names[x]
    y_label = axes_names[y]	
    

    """  Builds a list where each item is an array of a single GC """
    stuff = f.read_data(results_file)
    clus_data, data = [], []
    for i in range(len(stuff[:,0])):
        if (i==(len(stuff[:,0])-1) or (stuff[i,0] < stuff[(i-1),0]) ):
            #print sc.array(data), '\n'
            clus_data.append(sc.array(data))
            data = []
        row = stuff[i,:].tolist()
        #print row
        data.append(row)
    clus_data = clus_data[1:]  #gets rid of phantom first array

    """ Plotting Code """
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax2 = ax1.twiny()
    #print len(clus_data)
    #print clus_data
    for i in range(len(clus_data)):
        label = NAME[i]
        gc = clus_data[i]
        if (plot_errors==1):  y_errors = gc[:,e]
        else:  y_errors = None
        ax1.errorbar(gc[:,x], gc[:,y], yerr=y_errors, fmt=formats[i], ecolor='g',
                         marker='D', markersize=3)
        plt.text( (gc[0,x]+xoff), (gc[0,y]+yoff), label, fontsize=8 )
    #plot fit function
    #x_1 = sc.arange(5.0, 45.0, 0.5)
    #y_1 = sigmoid_function(x_1, [0.52058131, 11.99855162, 0.76338916])
    #plt.plot(x_1, y_1, 'k-.')
    #x_2 = sc.arange(5.0, 45.0, 0.5)
    #y_2 = linear_function(x_2, [0.0, 0.315])# [-0.00627752, 1.30478929])
    #y_2 = 0.36*sc.ones(len(x_2))
    #plt.plot(x_2, y_2, 'k-.')
    # Set up axes
    x_tick = ax1.get_xticks()
    ax1.set_xticks(x_tick)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    new_ticks = []
    for i in range(len(x_tick)):
        new_ticks.append(round( (5.0*m.log10(x_tick[i]/0.01) + 4.18), 1) )
    x_lim = ax1.get_xlim()
    new_lim = (5.0*m.log10(x_lim[0]/0.01) + 4.18), (5.0*m.log10(x_lim[1]/0.01) + 4.18)
    ax2 = ax1.twiny()
    ax2.set_xlim(new_lim)
    #ax2.set_xticks(new_ticks)
    plt.xticks(x_tick, new_ticks)
    plt.xticks(fontsize=10)
#    if (plot != 4):
#        plt.plot([np.ma.min(stuff[:,0]), np.ma.max(stuff[:,0])], [two_means[plot], two_means[plot]],
#            'k--')
    #now finds the average y of all points in dist. bins of size 2.0,
    # and outputs the set of points
    if (find_average == 1):
        avg_in = f.read_data(file_to_average)
        distances = avg_in[:,x]
        values = avg_in[:,y]
        margins = [9,11,13,15,17,19,21,23,25,27,29]
        points = []
        for i in range(len(margins)):
            holder = []
            for j in range(len(avg_in[:,0])):
                if ( (avg_in[j,0] > (margins[i]-2) ) and avg_in[j,0] < margins[i]):
                    holder.append(avg_in[j,y])
            data_point = [ (margins[i]-1.0), sc.mean(holder), sc.std(holder) ]
            points.append(data_point)
        averages = sc.array(points)
        plt.errorbar(averages[:,0], averages[:,1], yerr=averages[:,2], fmt='o', color='b',
                     ecolor='b', markersize=4)
        f.write_data(averages, averages_file,
                     header=('Distance, Average, Standard Deviation for '+axes_names[y]) )
        print '#-Printed and plotted averages and standard deviations'
    #plt.title('Plot of convolved data')
    ax1.set_xlabel(x_label)
    ax2.set_xlabel(r'$\bar{g}_0$')
    ax1.set_ylabel(y_label)
    fig.savefig(plot_file, papertype='letter')
    plt.close('all')
    print '#-plot successfully created'
print '#---All Plots sucessfully created'
