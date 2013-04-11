#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import gauss as ga

"""This program contains code for manipulating .csv files, 
esp. for Globular Cluster analysis.

Matthew Newby (RPI), Jan 14, 2009
"""

#These variable are to be changed based on the task:::
#"bool" imnlies a boolean value - 1 or 0,
GCname = "NGC 7089 con-23.2kpc f-turnoff"  #Name of cluster being analyzed
myfile = "convolved_7089_23.2kpc_cluster_f.csv"     #Name of input .csv file
distance = 11.5                #Distance to cluster in kpc
in_area = 0.0220        #Area of cluster data set 

background = 1        #subtract a sample background from the dataset - bool
backfile = "convolved_7089_23.2kpc_background_f.csv"    #background .csv file
back_area = 0.4120    #Area of background data set

write_txt_file = 0                 #create tab-delimited .txt file
txtfile = "test.txt"    #Name of .txt file to output, if any

cut_data_set = 0     #to cut the data set - bool
center_ra = 323.362    #ra value at center of cluster    
center_dec = -0.826    #dec value at center of cluster
inner_cut = 0.18    #inner radial boundry of dataset (degrees from center)
outer_cut = 10.0   #outer radial boundry of dataset (degrees from center) 
cut_txt_file = 0     #make a .txt file of the cut data - bool
cuttxtfile = "test_cut.txt"  #Name of cut data set to output, if any
cut_csv_file = 0     #make a .csv file of the cut data - bool
cutcsvfile = "test.csv"  #Name of cut data set to output, if any
plotcut = 0        #plot the cut data set on a scatter plot - bool

subplots = 0       #put all of the plots on one subplot -bool

set_bin_size = 1   #choose bin size manually, overrides nbins - bool
bin_sizes = 0.1    #set the size of each bin, in magnitudes
nbins = 100         #bins in histograms, below - might make this dynamic inside of code; 60 is good

cut_hist = 1       #cuts off bins outside of cut range, below. -bool
low_cut = 2.7    #low limit of histogram cut range - absolute mag.
high_cut = 8.7   #high limit of histogram cut range - absolute mag.
hist_low_cut = low_cut + 5.*(m.log10(distance*1000) - 1.)
hist_high_cut = high_cut + 5.*(m.log10(distance*1000) - 1.)

scatterplot = 0    #Ra - Dec scatterplot - bool
histogram = 0      #histogram, binned in g - bool
HRdiagram = 1      #H-R digram, g vs. g-r - bool
HRcontour = 0      #H-R diagram, cntour plot - bool

Nhistgauss = 1     #normalized histogram in g, with guassian fit overlay - bool
color_scatter = 1  #Creates a scatter plot using different colors for the background and cluster stars. -bool
plot_params = 1    #Includes parameters on the plot

HRplot_abs_mag = 1     #H-R plot uses absolute magnitude for y-axis instead of apparent - bool
set_HR_limits = 1      #H-R plot to use custom y-limits - bool
gridbins = 200         #H-R diagram density plot x bin number
HR_min = 7.0           #H-R custom minimum; use higher value for correct view
HR_max = 1.0          #H-R custom maximum; use lower value for correct view
HR_xmin = 0.0
HR_xmax = 0.6

convolve = 0           #Determines whether to convolve the data set or not; not used!

#Definitions:::

def run():
    i = 1  #figure iterator

    nrows = size(myfile)
    datarows = nrows - 1
     
    col1 = sc.zeros(datarows)  #RA
    col2 = sc.zeros(datarows)  #Dec
    col3 = sc.zeros(datarows)  #dered_g
    col4 = sc.zeros(datarows)  #dered_r

    col1, col2, col3, col4 = loadcsv(myfile, nrows)
    if background == 1:
        col5, col6, col7, col8 = loadcsv(backfile)
    
    if write_txt_file == 1:
        txt_success = writetxt(col1, col2, col3, col4, datarows, txtfile )
        if txt_success == 1:
            print 'original data successfully saved as', txtfile

#    if convolve == 1:
#        randomSeed = 138435
#        np.random.seed(randomSeed)
#        print "random seed:", randomSeed
#        con3 = convolve_column(col3, 'g')
#        con4 = convolve_column(col4, 'r')
#        col3 = con3
#        col4 = con4
#        print '-cluster g and r convolved'
        #if (background == 1):
            #con7 = convolve_column(col7, 'g')
            #col7 = con7
            #print '-background g convolved'

    if cut_data_set == 1:
        cola, colb, newsize = make_cut(col1, col2, col3, col4, datarows)
        if plotcut == 1:
            plt.figure(i)
            plt.scatter(cola, colb, marker='+')
            plt.axis('scaled')
            plt.title('Sky Plot of cut data from ' + GCname)
            plt.xlabel('RA')
            plt.ylabel('Dec')
            i = i + 1

    if subplots == 1:
        plt.figure(i)
        i = i + 1

    if scatterplot == 1:
        if subplots == 1:
            plt.subplot(221) 
        else:
            plt.figure(i)
            i = i + 1
        plt.scatter(col1, col2, marker='+')
        plt.axis('scaled')
        plt.title('Sky Plot of ' + GCname)
        plt.xlabel('RA')
        plt.ylabel('Dec')

    if histogram == 1:
        if subplots == 1:
            plt.subplot(222) 
        else:
            plt.figure(i)
            i = i + 1
 	numbins = 100
        plt.hist(col3, numbins, histtype='step')
        plt.title('Original Histogram of ' + GCname)
        plt.xlabel('g')
        plt.ylabel('counts')

    if Nhistgauss == 1:
        #x_hist, y_hist = normal_binned(col3, datarows)
        x_hist, y_hist = binned_background(col3, col7)
        nbins = len(x_hist)
        fit_params = ga.gauss_fit(x_hist, y_hist, nbins, 1,1)
        axis = 0.0 #(np.ma.max(x_hist) - np.ma.min(x_hist))/2.0 #non-zero lengthens x-axis
        x_line = sc.arange( (np.ma.min(x_hist)-axis), (np.ma.max(x_hist)+axis), 0.01)
        y_line = sc.zeros(len(x_line))
        mu_apparent = fit_params[0]
        fit_params[0] = fit_params[0] - 5.*(m.log10(distance*1000) - 1.)
        mu_str = '$\mu_M=$ ' + str(fit_params[0])
        sigmal_str = '$\sigma_l=$ ' + str(fit_params[1])
        sigmar_str = '$\sigma_r=$ ' + str(fit_params[2])
        A_str = '$A=$ ' + str(fit_params[3])
        mug_str = '$\mu_g=$ ' + str(mu_apparent)
        print '-best fit parameters:', fit_params
        print '-apparent magnitude mu:', mu_apparent
        fwhm = ga.full_width_half_max(fit_params[1], fit_params[2])
        FWHM_str = '$fwhm=$ ' + str(fwhm)
        print '-full width at half max: ', fwhm
        #Shifts x values to absolute magnitudes, creates a gaussian profile from params
        for j in range(len(x_hist)):
            x_hist[j] = x_hist[j] - 5.*(m.log10(distance*1000) - 1.)
        for j in range(len(x_line)):
            x_line[j] = x_line[j] - 5.*(m.log10(distance*1000) - 1.)
            if x_line[j] < fit_params[0]:
                sigma = fit_params[1]
            else:
                sigma = fit_params[2]
            exponent = -0.5 * ( (x_line[j] - fit_params[0]) / sigma )**2
            stuff = m.exp(exponent)
            y_line[j] = stuff*fit_params[3]
        #Creates a line to simulate binning
        x_bins = sc.arange( (np.ma.min(x_hist)-axis), (np.ma.max(x_hist)+axis), 0.01)
        y_bins = sc.zeros(len(x_bins))
        span = x_hist[1] - x_hist[0]
        bin_min = x_hist[0] + (span/2.0) #Starts at first bin-bin changover
        k=0
        for j in range(len(x_bins)):
            if (x_bins[j] > (bin_min + (span*k))):
                k = k + 1
            y_bins[j] = y_hist[k]

        #Shifts original x values to absolute magnitudes, prepares them for binning
        #x_shift = sc.zeros(datarows)
        #x_shift = col3 - 5.0*(m.log10(distance*1000) - 1.)
        if subplots == 1:
            plt.subplot(223) 
        else:
            plt.figure(i)
            i = i + 1
        plt.plot(x_bins, y_bins, 'b-')
        plt.scatter(x_hist, y_hist, marker='d' )
        #plt.hist(x_shift, nbins, histtype='step') 
        plt.plot(x_line, y_line, 'g-')
        plt.title('line fit of ' + GCname)
        plt.xlabel('Absolute Magnitude')
        plt.ylabel('Counts')
        if (plot_params == 1):
            txt_string = mu_str+'\n'+sigmal_str+'\n'+sigmar_str+'\n'+A_str+'\n'+mug_str+'\n'+FWHM_str
            plt.text(np.ma.min(x_line), (np.ma.max(y_line)*0.6), (txt_string))


    if HRdiagram == 1:
        gminusr = sc.zeros(datarows)
        gminusr = col3 - col4
        magnitudes = sc.zeros(datarows)
        if HRplot_abs_mag == 1:
            magnitudes = ( col3 - 5.*(m.log10(distance*1000) - 1.) )
        else:
            magnitudes = col3
        if subplots == 1:
            plt.subplot(224) 
        else:
            plt.figure(i)
            i = i + 1
        plt.scatter(gminusr, magnitudes, marker='+')
        if set_HR_limits == 1:
            plt.ylim(HR_min, HR_max)
            plt.xlim(HR_xmin, HR_xmax)
        else:
            M_min = np.ma.min(magnitudes)
            M_max = np.ma.max(magnitudes)
            plt.ylim(M_max, M_min)
        plt.title('HR diagram of ' + GCname)
        plt.xlabel('g-r')
        plt.ylabel('Mg')

    if HRcontour == 1:
        gminusr = sc.zeros(datarows)
        gminusr = col3 - col4
        magnitudes = sc.zeros(datarows)
        if HRplot_abs_mag == 1:
            magnitudes = ( col3 - 5.*(m.log10(distance*1000) - 1.) )
        else:
            magnitudes = col3
        plt.figure(i)
        i = i + 1
        M_min = np.ma.min(magnitudes)
        M_max = np.ma.max(magnitudes)
        plt.ylim(M_max, M_min)
        plt.hexbin(gminusr, magnitudes, gridsize=gridbins, bins='log', cmap=mpl.cm.gray_r)
        if set_HR_limits == 1:
            plt.ylim(HR_min, HR_max)
            plt.xlim(HR_xmin, HR_xmax)
        plt.title('HR diagram of ' + GCname)
        plt.xlabel('g-r')
        plt.ylabel('Mg')


    if color_scatter == 1:
        plt.figure(i)
        i = i + 1
        plt.scatter(col1, col2, c='blue', marker='+')
        plt.scatter(col5, col6, c='red', marker='x')
        plt.axis('scaled')
        plt.title('Sky Plot of ' + GCname)
        plt.xlabel('RA')
        plt.ylabel('Dec')


    plt.show()
    return '-Tasks successfully completed'


def size (infile):
    #opens a file and returns the number of rows in the file
    usef = open(infile, "rb")
    length = len( usef.readlines() )
    #see below for .csv version to find file length; not as nice as this one.
    print 'your file is', length, 'lines long'
    usef.close()
    return length


def loadcsv( infile, p=0, debug=0 ):
    #takes in a 4 column .csv file, with text first line (header) 
    #and subsequent lines are data,
    #and converts it to numerical arrays
    #p is the number of rows in the file
    #set debug to 1 to print out extra info

    #load file and read .csv format
    usefile = open(infile, "rb")

    #find the length of the file if no value for p given
    if p == 0:
        dataset = csv.reader(usefile, delimiter=',')
        for row in dataset:
            p = p + 1
        print '-Number of rows:', p
        usefile.seek(0)

    #Initialize arrays, points is the actual number of data points (0th line is header)
    points = p-1
    ra = sc.zeros(points)
    dec = sc.zeros(points)
    de_g = sc.zeros(points)
    de_r = sc.zeros(points)

    #resets file and rereads it (not sure if necessary)
    dataset = csv.reader(usefile, delimiter=',')

    #write file contents to memory
    i = 0
    for row in dataset:
        s_ra, s_dec, s_de_g, s_de_r = row
        if debug == 1:
            print row
        #pick off header 
        if i == 0:
            try:
                ra[i], dec[i], de_g[i], de_r[i] = eval(s_ra), eval(s_dec), eval(s_de_g), eval(s_de_r)
                print '-first line is not header'
                print '-please add a header line.  This program will now close'
                usefile.close()
            except NameError:
                hdr1, hdr2, hdr3, hdr4 = s_ra, s_dec, s_de_g, s_de_r
                print '-first line is header; data format:', hdr1, hdr2, hdr3, hdr4
            else:
                print '-something bad happened during file read...'
        #Write the data to all i (from 0 to points)
        else:
            ra[i-1], dec[i-1], de_g[i-1], de_r[i-1] = eval(s_ra), eval(s_dec), eval(s_de_g), eval(s_de_r)
            if debug == 1:
                print ra[i-1], dec[i-1], de_g[i-1], de_r[i-1]
        i = i + 1
        #Note that "i-1" is used so that the arrays fill from the 0th element up to "points"
   
    usefile.close()
    print '-Read', infile
    #returns the array values
    return ra, dec, de_g, de_r


def writecsv(ra, dec, dered_g, dered_r, points, filename):
    #write the dataset to a csv file
    fhandle = open(filename, 'w')
    outfile = csv.writer(fhandle, delimiter=',')
    for j in range(points):
        a, b, c, d = str(ra[j]), str(dec[j]), str(dered_g[j]), str(dered_r[j])
        e = [a, b, c, d]
        outfile.writerow( e )
    fhandle.close()
    return 1


def writetxt(ra, dec, dered_g, dered_r, points, filename ):
    #write a column-delimited ".txt" file with the data points in it
    fhandle = open(filename, 'w')
    for j in range(points):
        a, b, c, d = str(ra[j]), str(dec[j]), str(dered_g[j]), str(dered_r[j])
        fhandle.write( a + '\t' + b + '\t' + c + '\t' + d + '\n' )
    fhandle.close()
    return 1


def make_cut(col1, col2, col3, col4, points):
    #cuts a ring out of the data, cetered on the gc center 
    #save new data set to file(s), print out plot

    #find the data points that belong in the new set    
    holder = sc.zeros(points, int)
    in_set = 0
    for k in range(points):
        d = m.sqrt( (col1[k] - center_ra)**2 + (col2[k] - center_dec)**2 ) 
        if (d > inner_cut) and (d < outer_cut):
            holder[k] = 1
            in_set = in_set + 1
    outliers = (points - in_set)
    print 'Removed', outliers, 'from data set'

    #Make the new data set
    new1 = sc.zeros(in_set)
    new2 = sc.zeros(in_set)
    new3 = sc.zeros(in_set)
    new4 = sc.zeros(in_set)
    h = 0
    for l in range(points):
        if holder[l] == 1:
            new1[h], new2[h], new3[h], new4[h] = col1[l], col2[l], col3[l], col4[l]
            h = h + 1
    if h != in_set:
        print 'ERROR: not all cut data written properly!'

    #make files
    if cut_txt_file == 1:
        success1 = writetxt(new1, new2, new3, new4, in_set, cuttxtfile )
        if success1 == 1:
            print 'cut data successfully saved as', cuttxtfile

    if cut_csv_file == 1:
        success2 = writecsv(new1, new2, new3, new4, in_set, cutcsvfile )
        if success2 == 1:
            print 'cut data successfully saved as', cutcsvfile


    return new1, new2, in_set


def binned_background(data_x, back_x):
    #bins the data and the background file, then subtracts 
    #the background from the data

    #Makes the range go from the minimum of both data set to the maximum
    data_min = np.ma.min(data_x)
    back_min = np.ma.min(back_x)
    if (data_min <= back_min):
        x_min = data_min
    else:
        x_min = back_min
    data_max = np.ma.max(data_x)
    back_max = np.ma.max(back_x)
    if (data_max >= back_max):
        x_max = data_max
    else:
        x_max = back_max
    if set_bin_size == 1:
        bin_size = bin_sizes
        nbins = int( ((x_max+bin_size) - x_min)/bin_size )  #add a bin to compensate for rounding
    else:
        bin_size = ( (x_max - x_min)/ nbins )
    print '-', nbins, 'bins of size', bin_size, 'magnitudes'

    #Creates the bins
    bin_centers = sc.zeros(nbins)
    bin_counts = sc.zeros(nbins)
    for j in range(nbins):
        bin_centers[j] = x_min + (j*bin_size) + (bin_size*0.5)

    #Bin data, and subtract background; no bin counts below zero
    u = 0  #Tracks the number of binned points from dataset
    v = 0  #Tracks the number of "unbinned" points, taken from background
    w = 0.0  #Tracks number of points removed due to hist_cut selection
    #bins all values between bin max and bin min, including points at bin min
    for k in range(nbins):
        bin_min = x_min + (k*bin_size)
        bin_max = x_min + ((k+1)*bin_size)
        for l in range(len(data_x)):
            if (bin_min <= data_x[l] < bin_max):
                bin_counts[k] = bin_counts[k] + 1.0
                u = u + 1  #bug test 
    #Bins the last data point
    for k in range(len(data_x)):
        if (data_x[k] == x_max):  
            bin_counts[-1] = bin_counts[-1] + 1.0
            u = u + 1  #tracking
    #Removes background from bins
    for k in range(nbins):
        bin_min = x_min + (k*bin_size)
        bin_max = x_min + ((k+1)*bin_size)
        back_bin = 0.0
        for l in range(len(back_x)):
            if (bin_min <= back_x[l] < bin_max):
                back_bin = back_bin + 1.0
                v = v + 1  #tracking
        #Compensating for area differences
        back_bin = back_bin*(in_area/back_area)
        bin_counts[k] = bin_counts[k] - back_bin
        if (bin_counts[k] < 0):
            bin_counts[k] = 0.000
    #Removes the last data point, if necessary
    for k in range(len(back_x)):
        if (back_x[k] == x_max):  
            bin_counts[-1] = bin_counts[-1] - (in_area/back_area)
            v = v + 1  #tracking
        if (bin_counts[-1] < 0):
            bin_counts[-1] = 0.00
    #Cuts the data set, if selected
    if (cut_hist == 1):
        print '-Using only data between M=', low_cut, 'and', high_cut
        for k in range(nbins):
            if (bin_centers[k] < hist_low_cut):
                w = w + bin_counts[k]  #tracking
                bin_counts[k] = 0.0
            if (bin_centers[k] > hist_high_cut):
                w = w + bin_counts[k]  #tracking
                bin_counts[k] = 0.0

    print '-Total Data in:', len(data_x), '; total points binned:', u
    print '-Total Background in:', len(back_x), '; total background removed:', v
    print '-Total number of points removed due to hist cuts:', w
    return bin_centers, bin_counts

def normal_binned(x_val, points):
    #x values are dered_g; points is the number of data points in x
    #Creates an x vs y data set, with x in bins and y in counts per bin.

    g_min = np.ma.min(x_val)
    g_max = np.ma.max(x_val)
    bin_size = ( (g_max - g_min)/(nbins) )

    bin_centers = sc.zeros(nbins)
    bin_counts = sc.zeros(nbins, int)
    for j in range(nbins):
        bin_centers[j] = g_min + (j*bin_size) + (bin_size*0.5)

    z = 0 #bug test
    for k in range(nbins):
        bin_min = g_min + (k*bin_size)
        bin_max = g_min + ((k+1)*bin_size)
        for l in range(points):
            if (bin_min <= x_val[l] < bin_max):
                bin_counts[k] = bin_counts[k] + 1
                z = z + 1  #bug test 
    for n in range(points):
        if (x_val[n] == g_max):  
            bin_counts[-1] = bin_counts[-1] + 1
            z = z + 1  #bug test 

    #debugging crap
    #for j in range(nbins):
    #    print bin_centers[j], bin_counts[j]
    print '-total data in:', points, '; total points binned:', z

    #bin cut
   
    return bin_centers, bin_counts


def R_squared_gauss (mu, sigma_left, sigma_right, data_x, data_y, points):
    #This function returns the r-squared value of a data set, x,y, that is "points" long,
    #and the current gaussian parameter set, mu, sigma.  Note that this gaussian is double-sided.

    line_y = sc.zeros(points)
    
    #Create the gaussian from parameter sets
    for j in range(points):
        if data_x[j] < mu:
            sigma = sigma_left
        else:
            sigma = sigma_right
        normalizer = (np.sqrt(2*m.pi))*sigma
        exponent = -0.5 * ( (data_x[j] - mu) / sigma )**2
        stuff = m.exp(exponent)
        line_y[j] = stuff / normalizer

    #Find r-squared value
    r_squared = 0.0
    for k in range(points):
        diff = (data_y[k] - line_y[k])**2
        r_squared = r_squared + diff

    return r_squared

def convolve_column(column, color):
    a_g = 0.0
    b_g = 0.790391
    c_g = -19.8928
    a_r = 0.0
    b_r = 0.766309
    c_r = -19.0334
    near_d = 7.5
    far_d = 23.2
    diff_mag = 5.0*(np.log10(far_d / near_d))
    l = len(column)
    new_column = sc.zeros(l)
    for i in range(l):
        if color == 'g':
            sigma_con = a_g + np.exp(b_g*(column[i] + diff_mag) + c_g)
        elif color == 'r':
            sigma_con = a_r + np.exp(b_r*(column[i] + diff_mag) + c_r)
        else:  
            sigma_con = 0.0
            print '!!!   Invalid Color Selection   !!!'
        new_column[i] = np.random.normal(column[i], sigma_con)
    return new_column

