import math as m
import numpy as np
import scipy as sc
import files as f
import gctools_enc as gce
import gauss_fit_enc as gfe
import GCplot_enc as gcp

'''python script for running convolve analysis scripts
Matthew Newby, June 6, 2010'''

datafile = 'HR_Pal_5_cluster.csv'
NAME = ['Pal5_ah23', 'Pal5_ah24', 'Pal5_ah25', 'Pal5_ah26', 'Pal5_ah27', 'Pal5_ah28']
CONVOLVE = [23.2, 24.2, 25.2, 26.2, 27.2, 28.2]
real_dist = 23.2
repeat_runs = 3

print '#', CONVOLVE
print '#-Total iterations:', repeat_runs
for i in range(repeat_runs):
    for run in range(len(NAME)):
        gc_name = NAME[run]
        con_dist = CONVOLVE[run]
        print '#*****  Starting run', gc_name, 'from distace', real_dist, \
        'to distance:', con_dist ,'*****'
        main_data = f.read_csv(datafile)
        con_clus = gce.convolve(main_data, real_dist, con_dist)
        if gce.con_analysis(main_data, con_clus, right_limit=0.3, left_limit=0.1) == 1:
            print '#*****', gc_name, '-process complete*****'
print '#---ALL RUNS COMPLETE---'
