import math as m
import numpy as np
import scipy as sc
import files as f
import gctools_enc as gce
import time

'''python script for running gc scripts
Matthew Newby, Feb 16, 2011'''

def dummy ():
    return 0

NAME = ['4147', '5024', '5053', '5272', '5466', '5904', '6205', '6341', '7078', '7089', 'Pal5']
SKIP = [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1]  #Set to '1' to skip a cluster
STEP = 5.0 
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
DISTANCE = [19.3, 18.7, 18.5, 10.4, 15.6, 8.0, 7.7, 8.7, 11.0, 11.5, 21.0] 
AREA_IN = [0.0145, 0.1679, 0.0537, 0.1788, 0.0826, 0.076, 0.0664, 0.0175, 0.0504, 0.022, 0.0238]
AREA_BACK = [0.0915, 1.1548, 0.6922, 2.1823, 0.9737, 1.8537, 3.2279, 0.0635, 0.2517, 0.4120, 0.2740]
BIN_SIZE = 0.2
THRESHOLD = None #10
cut_hist = 1
low_cut = 2.0 #2.5
high_cut = 8.0
HIGH_CUTS = []  #use this one when cuts are the same for all
#HIGH_CUTS = [6.0, 6.0, 6.0, 6.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]  #5272, short series
#HIGH_CUTS = [6.3, 6.3, 6.5, 6.5, 7.1, 8.0, 8.0, 8.0, 8.0, 8.0]  #7089, short series
AVERAGES = 10
BLUE_LIMIT = 0.1
RED_LIMIT = 0.3
plot_file = 1
np.random.seed(int(time.time()))
print '#-Random seed =', int(time.time())

def convolve(data_in, real_dist, con_dist):
    a_g = 0.0
    b_g = 0.790391
    c_g = -19.8928
    a_r = 0.0
    b_r = 0.766309
    c_r = -19.0334
    diff_mag = 5.0*(np.log10(con_dist / real_dist))
#    print '#---from ', real_dist, 'kpc to', con_dist, 'kpc'
    length = len(data_in)
    g = sc.zeros(length)
    r = sc.zeros(length)
    for i in range(length):
        sigma_g = np.sqrt((np.exp(b_g*(data_in[i][2]+diff_mag))**2 - np.exp(b_g*data_in[i][2])**2))*np.exp(c_g)
        print 'g', sigma_g
        sigma_r = np.sqrt((np.exp(b_r*(data_in[i][3]+diff_mag))**2 - np.exp(b_r*data_in[i][3])**2))*np.exp(c_r)
        print 'r', sigma_r
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

data_out = []
for cluster in range(len(NAME)):
    if (SKIP[cluster] == 1):  continue
    con_dist = DISTANCE[cluster]
    real_dist = DISTANCE[cluster]
    gc_name = NAME[cluster]
    main_data = f.read_csv(datafile[cluster])
    l, w = main_data.shape
    print '#-Starting run', gc_name
    print '#-Total number of stars in cluster:', l
    ''' Sort the Data '''
    blue, yellow, red = [], [], []
    for i in range(l):
        """Remove stars from outside interesting range"""
        if (  (main_data[i,2] - 5.*(m.log10(real_dist*1000) - 1.) ) < low_cut ): continue
        if (  (main_data[i,2] - 5.*(m.log10(real_dist*1000) - 1.) ) > high_cut ): continue
        g_minus_r = (main_data[i,2] - main_data[i,3])
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
        '''Average all the values in 'values' '''
        hold_array = sc.array(values)
        holder = []
        for i in range(10):
            holder.append(sc.mean(hold_array[:,i]))
        print holder
        data_out.append(holder)
f.write_data(sc.array(data_out), 'new_output.txt')
print '#---All Done!'