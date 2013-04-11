import math as m
import numpy as np
import scipy as sc
import files as f
import gctools_enc as gce
import time

'''python script for running gc scripts
Matthew Newby, Feb 16, 2011'''

NAME = ['6205']
SKIP = [0]  #Set to '1' to skip a cluster
STEP = 1.0 
datafile = ['NGC_6205_all_cluster.csv']
#backfile = ['NGC_6205_all_background.csv'] #Not used
DISTANCE = [1.0]  #Psuedo-real distance; i.e., starting distance for convolutions
ACTUAL = [7.7]
THRESHOLD = None #10
cut_hist = 1
low_cut = 2.0
high_cut = 8.0
HIGH_CUTS = []  #use this one when cuts are the same for all
#HIGH_CUTS = [6.0, 6.0, 6.0, 6.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]  #5272, short series
#HIGH_CUTS = [6.3, 6.3, 6.5, 6.5, 7.1, 8.0, 8.0, 8.0, 8.0, 8.0]  #7089, short series
AVERAGES = 100 #00
BLUE_LIMIT = 0.2 #0.1
RED_LIMIT = 1.0 #0.3
plot_file = 1
ug_cut = 1  #1 to enforce u-g limit, 0 otherwise
det_eff = 0 # 1 to enforce detection efficiency, 0 otherwise.
np.random.seed(int(time.time()))
print '#-Random seed =', int(time.time())

def sigmoid_error(x, modulus=None):
    s = [0.9402, 1.6171, 23.5877]
    if modulus != None:
        s[2] = s[2] - modulus
    detection_efficiency = s[0] / (np.exp(s[1]*(x - s[2])) + 1.)
    return detection_efficiency

def convolve(data_in, real_dist, con_dist):
    """ shifts the u,g,r values of a data set by random normal values appropriate
        to the errors.  Stars remain at roughly the same 'distance,' but have
        their colors shifted"""
    a_u = 2.71e-03
    b_u = 0.80
    c_u = -19.2
    a_g = 3.11e-04
    b_g = 0.79
    c_g = -20.0
    a_r = -2.63e-05
    b_r = 0.80
    c_r = -19.8
    diff_mag = 5.0*(np.log10(con_dist / real_dist))
    #print '#-Convolving data with random seed =', randomSeed
    #print '#---from ', real_dist, 'kpc to', con_dist, 'kpc'
    length = len(data_in)
    u = sc.zeros(length)
    g = sc.zeros(length)
    r = sc.zeros(length)
    for i in range(length):
        #sigma_g = np.sqrt((np.exp(b_g*(data_in[i,2]+diff_mag))**2 - np.exp(b_g*data_in[i,2])**2))*np.exp(c_g)
        #sigma_r = np.sqrt((np.exp(b_r*(data_in[i,3]+diff_mag))**2 - np.exp(b_r*data_in[i,3])**2))*np.exp(c_r)
        sigma_u = np.sqrt(( (a_u + np.exp(b_u*(data_in[i][4]+diff_mag)+c_u))**2 - (a_u + np.exp((b_u*data_in[i][4]) + c_u))**2))
        sigma_g = np.sqrt(( (a_g + np.exp(b_g*(data_in[i][2]+diff_mag)+c_g))**2 - (a_g + np.exp((b_g*data_in[i][2]) + c_g))**2))
        sigma_r = np.sqrt(( (a_r + np.exp(b_r*(data_in[i][3]+diff_mag)+c_r))**2 - (a_r + np.exp((b_r*data_in[i][3]) + c_r))**2))
        if (sigma_g > 0.0):
            g[i] = np.random.normal(data_in[i][2], sigma_g)
        else:  
            g[i] = data_in[i][2]
        if (sigma_r > 0.0):
            r[i] = np.random.normal(data_in[i][3], sigma_r)
        else:  
            r[i] = data_in[i][3]
        if (sigma_u > 0.0):
            u[i] = np.random.normal(data_in[i][4], sigma_u)
        else:  
            u[i] = data_in[i][4]
    m = sc.zeros((length,5))
    for i in range(length):
        m[i,0] = data_in[i][0]
        m[i,1] = data_in[i][1]
        m[i,2] = g[i]
        m[i,3] = r[i]
        m[i,4] = u[i]
    return m

""" MAIN """
if __name__ == "__main__":
    data_out = []
    for cluster in range(len(NAME)):
        if (SKIP[0] == 1):  continue
        con_dist = DISTANCE[0]
        real_dist = DISTANCE[0] #pseudo-real distance; starting distance
        actual_dist = ACTUAL[0]
        gc_name = NAME[0]
        main_data = f.read_csv(datafile[0])
        l, w = main_data.shape
        """ Moduli for real cluster and con cluster """
        mod_actual = 5.*(m.log10(actual_dist*1000) - 1.)
        mod_real = 5.*(m.log10(real_dist*1000) - 1.)
        """ Make the data consistant with observation at starting distance """
        main_data[:,2] = main_data[:,2] - (mod_actual - mod_real) #- 4.515
        main_data[:,3] = main_data[:,3] - (mod_actual - mod_real) #- 4.515
        print '#-Starting run', gc_name
        print '#-Total number of stars in cluster:', l
        ''' Sort the Data '''
        blue, yellow, red = [], [], []
        for i in range(l):
            """Remove stars from outside interesting range"""
            if (  (main_data[i,2] - 5.*(m.log10(real_dist*1000) - 1.) ) < low_cut ): continue
            if (  (main_data[i,2] - 5.*(m.log10(real_dist*1000) - 1.) ) > high_cut ): continue
            """Sort stars accordingly """
            g_minus_r = (main_data[i,2] - main_data[i,3])
            if ug_cut == 1:
                if (main_data[i,4] - main_data[i,2]) < 0.4:continue
            if g_minus_r < BLUE_LIMIT:
                blue.append(main_data[i,:])
            elif g_minus_r > RED_LIMIT:
                red.append(main_data[i,:])
            else:
                yellow.append(main_data[i,:])
        base_blue, base_red, base_yellow = len(blue), len(red), len(yellow)
        #Output array is rows of data at each d_{eff}; 
        data_out.append([real_dist, base_blue, 0, 0, 0, base_yellow, 0, 0, 0, base_red])
        print '#-Original number of stars in blue, yellow, and red regions:', base_blue, base_yellow, base_red
        while (con_dist < 80.0):
            con_dist = con_dist + STEP
            values = []
            for iteration in range(AVERAGES):
                killed = 0
                holder = []
                holder.append(con_dist)
                '''CONVOLVE DATA'''
                print '#---convolving from', real_dist, 'kpc to', con_dist, 'kpc'
                blue_con = convolve(blue, real_dist, con_dist)
                red_con = convolve(red, real_dist, con_dist)
                yellow_con = convolve(yellow, real_dist, con_dist)
                '''Count Data - blue'''
                blue_stars, yellow_stars, red_stars = 0, 0, 0
                for i in range(len(blue_con)):
                    g_minus_r = (blue_con[i][2] - blue_con[i][3])
                    if ug_cut == 1:
                        if (blue_con[i][4] - blue_con[i][2]) < 0.4:  continue
                    if det_eff == 1:
                        if sigmoid_error(blue_con[i][2]) < np.random.uniform():
                            killed = killed + 1
                            continue
                    if g_minus_r < BLUE_LIMIT:
                        blue_stars = blue_stars + 1
                    elif g_minus_r > RED_LIMIT:
                        red_stars = red_stars + 1
                    else:
                        yellow_stars = yellow_stars + 1
                holder.append(blue_stars)
                holder.append(yellow_stars)
                holder.append(red_stars)
                '''Count Data - yellow'''
                blue_stars, yellow_stars, red_stars = 0, 0, 0
                for i in range(len(yellow_con)):
                    g_minus_r = (yellow_con[i][2] - yellow_con[i][3])
                    if ug_cut == 1:
                        if (yellow_con[i][4] - yellow_con[i][2]) < 0.4:  continue
                    if det_eff == 1:
                        if sigmoid_error(yellow_con[i][2]) < np.random.uniform():
                            killed = killed + 1
                            continue
                    if g_minus_r < BLUE_LIMIT:
                        blue_stars = blue_stars + 1
                    elif g_minus_r > RED_LIMIT:
                        red_stars = red_stars + 1
                    else:
                        yellow_stars = yellow_stars + 1
                holder.append(blue_stars)
                holder.append(yellow_stars)
                holder.append(red_stars)
                '''Count Data - red'''
                blue_stars, yellow_stars, red_stars = 0, 0, 0
                for i in range(len(red_con)):
                    g_minus_r = (red_con[i][2] - red_con[i][3])
                    if ug_cut == 1:
                        if (red_con[i][4] - red_con[i][2]) < 0.4:  continue
                    if det_eff == 1:
                        if sigmoid_error(red_con[i][2]) < np.random.uniform():
                            killed = killed + 1
                            continue
                    if g_minus_r < BLUE_LIMIT:
                        blue_stars = blue_stars + 1
                    elif g_minus_r > RED_LIMIT:
                        red_stars = red_stars + 1
                    else:
                        yellow_stars = yellow_stars + 1
                holder.append(blue_stars)
                holder.append(yellow_stars)
                holder.append(red_stars)
                '''stick it in the value holder'''
                print '#---', holder
                values.append(holder)
                print "### DEAD", killed, "###"
            '''Average all the values in 'values' '''
            hold_array = sc.array(values)
            holder = []
            for i in range(10):
                holder.append(sc.mean(hold_array[:,i]))
            print holder
            data_out.append(holder)
    f.write_data(sc.array(data_out), 'new_output.txt', header="distance (kpc), blue-in-blue, b-y, b-r; y-b, y-y, y-r; r-b, r-y, r-r")
    print '#---All Done!'

"""  Probably crap!
def convolve(data_in, real_dist, con_dist):
    a_u = 2.71e-03
    b_u = 0.80
    c_u = -19.2
    a_g = 3.11e-04
    b_g = 0.79
    c_g = -20.0
    a_r = -2.63e-05
    b_r = 0.80
    c_r = -19.8
    diff_mag = 5.0*(np.log10(con_dist / real_dist))
#    print '#---from ', real_dist, 'kpc to', con_dist, 'kpc'
    length = len(data_in)
    g = sc.zeros(length)
    r = sc.zeros(length)
    u = sc.zeros(length)
    for i in range(length):
        sigma_g = np.sqrt((np.exp(b_g*(data_in[i][2]+diff_mag))**2 - np.exp(b_g*data_in[i][2])**2))*np.exp(c_g)
        #print 'g', sigma_g
        sigma_r = np.sqrt((np.exp(b_r*(data_in[i][3]+diff_mag))**2 - np.exp(b_r*data_in[i][3])**2))*np.exp(c_r)
        #print 'r', sigma_r
        sigma_u = np.sqrt(( (a_u + np.exp(b_u*(data_in[i,4]+diff_mag)+c_u))**2 - (a_u + np.exp((b_u*data_in[i,4]) + c_u))**2))
        if (sigma_g > 0.0):
            g[i] = np.random.normal(data_in[i][2], sigma_g)
        else:  
            g[i] = data_in[i][2]
        if (sigma_r > 0.0):
            r[i] = np.random.normal(data_in[i][3], sigma_r)
        else:  
            r[i] = data_in[i][3]
    m = sc.zeros((length,4))
    for i in range(length):
        m[i,0] = data_in[i][0]
        m[i,1] = data_in[i][1]
        m[i,2] = g[i]
        m[i,3] = r[i]
    return m
"""