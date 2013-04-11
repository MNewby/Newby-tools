import math as m
import numpy as np
import scipy as sc
import files as f
import gctools_enc as gce
import gauss_fit_enc as gfe
import GCplot_enc as gcp

'''python script for running gc scripts
Matthew Newby, June 6, 2010'''

NAME = ['5904_tO','5904_t9', '5904_t11','5904_t13','5904_t15','5904_t17','5904_t19','5904_t21','5904_t23']
CONVOLVE = [7.3, 9.3, 11.3, 13.3, 15.3, 17.3, 19.3, 21.3, 23.2]
datafile = 'HR_NGC_5904_cluster.csv'
backfile = 'HR_NGC_5904_background.csv'
real_dist = 7.3
in_area = 0.076
back_area = 1.8537
cut_hist = 1
low_cut = 2.0
high_cut = 7.0
plot_file = 1

for run in range(len(NAME)):
    gc_name = NAME[run]
    con_dist = CONVOLVE[run]
    print '*****  Starting run', gc_name, '*****'

    main_data = f.read_csv(datafile)
    back_data = f.read_csv(backfile)
    con_clus = gce.convolve(main_data, real_dist, con_dist)
    con_back = gce.convolve(back_data, real_dist, con_dist)
    con_clus_f = gce.select_stars(con_clus)
    con_back_f = gce.select_stars(con_back)

    lc,wc = con_clus_f.shape
    lb,wb = con_back_f.shape
    cluster_x = sc.zeros(lc)
    background_x = sc.zeros(lb)
    for i in range(lc):
        cluster_x[i] = ( con_clus_f[i,2] - 5.*(m.log10(real_dist*1000) - 1.) )
    for i in range(lb):
        background_x[i] = ( con_back_f[i,2] - 5.*(m.log10(real_dist*1000) - 1.) )

    hist_data = gce.binned_background(cluster_x, background_x, 0.1, in_area, back_area, 
        cut_hist, low_cut, high_cut)

    start_params = gfe.make_start_params(hist_data)

    best_RR1, best_parameters1 = gfe.do_MCMC(hist_data, start_params)
    print '-MCMC1 best:', best_RR1, best_parameters1
    print '-MCMC1 app mu:', (best_parameters1[0] + 5.*(m.log10(real_dist*1000) - 1.))
    print '-MCMC1 fwhm:', (gfe.full_width_half_max(best_parameters1[1],best_parameters1[2]))

    best_RR2, best_parameters2 = gfe.do_MCMC(hist_data, start_params)
    print '-MCMC2 best:', best_RR2, best_parameters2
    print '-MCMC2 app mu:', (best_parameters2[0] + 5.*(m.log10(real_dist*1000) - 1.))
    print '-MCMC2 fwhm:', (gfe.full_width_half_max(best_parameters2[1],best_parameters2[2]))

    if (best_RR1 < best_RR2):
        gd_start = best_parameters1
    else:
        gd_start = best_parameters2

    best_parameters = gfe.Gradient_Descent(hist_data, gd_start)
    best_RR = gfe.R_squared_gauss(hist_data, best_parameters)
    print '-GD best:', best_RR, best_parameters
    print '-GD app mu:', (best_parameters[0] + 5.*(m.log10(real_dist*1000) - 1.))
    print '-GD fwhm:', (gfe.full_width_half_max(best_parameters[1],best_parameters[2]))

    if (gcp.plot_hist(hist_data, gc_name, 1, best_parameters, 1, plot_file) == 1):
        print '-Histogram plotted'
    if (gcp.plot_HR_diagram(con_clus, real_dist, gc_name,[1.0,7.0],[0.0,0.6], plot_file) == 1):
        print '-HR diagram plotted'
    if (gcp.plot_infiles(main_data, back_data, gc_name, plot_file) == 1):
        print '-Skyplot finished'

    print '******', gc_name, '-process complete*****'
print '---ALL RUNS COMPLETE---'
