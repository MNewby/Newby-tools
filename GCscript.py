import math as m
import numpy as np
import scipy as sc
import files as f
import gctools_enc as gce
import gauss_fit_enc as gfe
import GCplot_enc as gcp
import Hessian_errors as hes

'''python script for running gc scripts
Matthew Newby, June 6, 2010'''

def dummy ():
    return 0 

NAME = ['4147', '5024', '5053', '5272', '5466', '5904', '6205', '6341', '7078', '7089', 'Pal5']
SKIP = [1,1,1,0,1,1,1,1,1,1,1]  #Set to '1' to skip a cluster
CONVOLVE = [ #[19.3], [18.7], [18.5], [10.4, 44.0], [15.6], [8.0], [7.7], [8.7], [11.0], [11.5], [21.0]] #singles 
#[19.3, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0],
[19.3, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
#[18.7, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0],
[18.7, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
#[18.5, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0],
[18.5, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
#[10.4, 15.4, 20.4, 25.4, 30.4, 35.4, 40.4, 45.4, 50.4, 55.4, 60.4],
[10.4, 44.0], #[10.4, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
[15.6, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
#[8.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0],
[8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
#[7.7, 10.7, 15.7, 20.7, 25.7, 30.7, 35.7, 40.7, 45.7, 50.7, 55.7, 60.7],
[7.7, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
[8.7, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
[11.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
[11.5, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0],
[21.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0]
] 
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
#HIGH_CUTS = [6.0, 6.0, 6.0, 6.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]  #5272, short series
#HIGH_CUTS = [6.3, 6.3, 6.5, 6.5, 7.1, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]  #7089, short series
plot_file = 1


for cluster in range(len(NAME)):
    if (SKIP[cluster] == 1):  continue
    for run in range(len(CONVOLVE[cluster])):
   
        gc_name = NAME[cluster] +'_'+ str(CONVOLVE[cluster][run])
        con_dist = CONVOLVE[cluster][run]
        real_dist = DISTANCE[cluster]
        in_area = AREA_IN[cluster]
        back_area = AREA_BACK[cluster]
        if (HIGH_CUTS != []):  high_cut = HIGH_CUTS[run]
        print '*****  Starting run', gc_name, '*****'

        main_data = f.read_csv(datafile[cluster])
        back_data = f.read_csv(backfile[cluster])

        con_clus = gce.convolve(main_data, real_dist, con_dist)
        #don't convolve background!  Background is static!!!
        #Actually, DO convolve background!  need to accurately remove it from cluster!
        # If cluster moves before subtraction, background needs to move, too!
        con_back = gce.convolve(back_data, real_dist, con_dist)  

        con_clus_f = gce.select_stars(con_clus, low_limit = 0.1)
        con_back_f = gce.select_stars(back_data, low_limit = 0.1)  

        lc,wc = con_clus_f.shape
        lb,wb = con_back_f.shape
        cluster_x = sc.zeros(lc)
        background_x = sc.zeros(lb)
        #modulus = None
        modulus = 5.*(m.log10(real_dist*1000) - 1.)
        for i in range(lc):
            cluster_x[i] = ( con_clus_f[i,2] - modulus )
        for i in range(lb):
            background_x[i] = ( con_back_f[i,2] - modulus )

        #Bin cluster and background and subtract; now also compensates for initial detection errors.
        hist_data = gce.binned_background(cluster_x, background_x, BIN_SIZE, in_area,
                                          back_area, cut_hist, low_cut, high_cut,
                                          THRESHOLD, modulus)
        
        far_modulus = 5.*(m.log10(con_dist*1000) - 1.)
        for bin in range(len(hist_data[:,0])):
            '''For Cluster Distributions'''
            hist_data[bin,1] = hist_data[bin,1] * gce.parabola_error(hist_data[bin,0], far_modulus)
            print hist_data[bin,1], gce.parabola_error(hist_data[bin,0])
            '''For Background distributions '''
            #hist_data[bin,1] = hist_data[bin,1] * gce.sigmoid_error(hist_data[bin,0], far_modulus)
            #print hist_data[bin,1], gce.sigmoid_error(hist_data[bin,0])
        #hist_data = gce.sigmoid_errors(hist_raw, real_dist, con_dist)
        """ """
        start_params = gfe.make_start_params(hist_data)

        best_RR1, best_parameters1 = gfe.do_MCMC(hist_data, start_params)
        print '#-MCMC1 best:', best_RR1, best_parameters1
        print '#-MCMC1 app mu:', (best_parameters1[0] + 5.*(m.log10(real_dist*1000) - 1.))
        print '#-MCMC1 fwhm:', (gfe.full_width_half_max(best_parameters1[1],best_parameters1[2]))

        best_RR2, best_parameters2 = gfe.do_MCMC(hist_data, start_params)
        print '#-MCMC2 best:', best_RR2, best_parameters2
        print '#-MCMC2 app mu:', (best_parameters2[0] + 5.*(m.log10(real_dist*1000) - 1.))
        print '#-MCMC2 fwhm:', (gfe.full_width_half_max(best_parameters2[1],best_parameters2[2]))

        if (best_RR1 < best_RR2):
            GDstart = best_parameters1
        else:
            GDstart = best_parameters2

        best_parameters = gfe.Gradient_Descent(hist_data, GDstart)
        best_RR = gfe.R_squared_gauss(hist_data, best_parameters)
        print '#-GD best:', best_RR, best_parameters
        print '#-GD app mu:', (best_parameters[0] + 5.*(m.log10(real_dist*1000) - 1.))
        print '#-GD fwhm:', (gfe.full_width_half_max(best_parameters[1],best_parameters[2]))

        print '#-GD best:', best_RR, best_parameters
        print '#-GD app mu:', (best_parameters[0] + 5.*(m.log10(real_dist*1000) - 1.))
        print '#-GD fwhm:', (gfe.full_width_half_max(best_parameters[1],best_parameters[2]))

        lh,wh = hist_data.shape
        sig_in = sc.ones(lh)
        step_size = [0.01, 0.001, 0.001, 0.01]
        errors = hes.get_hessian_errors(gfe.get_2gauss_y, best_parameters, hist_data[:,0],
                                         hist_data[:,1], sig_in, step_size)

        print '#- Distance, fwhm, mu, sigma_l, sigma_r, amplitude, errors for run', gc_name
        print '$ ', con_dist, (gfe.full_width_half_max(best_parameters[1],best_parameters[2])), \
            best_parameters[0], best_parameters[1], best_parameters[2], best_parameters[3], \
            errors[0], errors[1], errors[2], errors[3]
        
        if (gcp.plot_hist(hist_data, gc_name, 1, best_parameters, 0, plot_file) == 1):  
            print '#-Histogram plotted'
        if (gcp.plot_HR_diagram(con_clus, real_dist, gc_name,[1.0,7.0],[0.0,0.6], plot_file) == 1):
        #if (gcp.plot_HR_diagram(gce.select_stars_with_sigmoid(con_clus, con_dist, low_limit=0.0,high_limit=0.6),
        #                        real_dist, gc_name,[1.0,7.0],[0.0,0.6], plot_file) == 1):
            print '#-HR diagram plotted'
        if run == 0:
            if (gcp.plot_infiles(main_data, back_data, gc_name, plot_file) == 1):
                print '-Skyplot finished'
        print '******', gc_name, '-process complete*****'
        
        """ EXTRA"""
        #print cluster, run
        f.write_data(hist_data, (NAME[cluster]+'_'+str(round(CONVOLVE[cluster][run]))+'_'+'.txt'),
                     delimiter='\t', header=str(best_parameters))
print '#---ALL RUNS COMPLETE---'
