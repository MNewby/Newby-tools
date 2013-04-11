import math as m
import numpy as np
import scipy as sc
import files as f

'''python script for quickly averaging data.
Matthew Newby, November 15,2010'''

all_data = f.read_data('best_results.txt')
point_data = f.read_data('best_results_firstpoint.txt')
#get mu averages:
print '#-average in mu of first points:', sc.mean(point_data[:,2]), sc.std(point_data[:,2])
print '#-average in mu of all points:', sc.mean(all_data[:,2]), sc.std(all_data[:,2])
#get sig_l averages:
print '#-average in sig_l of first points:', sc.mean(point_data[:,3]), sc.std(point_data[:,3])
print '#-average in sig_l of all points:', sc.mean(all_data[:,3]), sc.std(all_data[:,3])