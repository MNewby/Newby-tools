import scipy as sc
import math as m
import numpy as np
import files as f
import GCplot_enc as gcp

"""Script for quick plotting of HR diagrams.

Matthew Newby, RPI, Aug 29, 2010"""

gc_name = 'NGC_6341'
clus_file = 'HR_' + gc_name + '_background.csv'
back_file = 'HR_' + gc_name + '_cluster.csv'
#load file
clus_data = f.read_csv(clus_file)
back_data = f.read_csv(back_file)
#plotfile
if (gcp.plot_infiles(back_data, clus_data, gc_name)==1):
    print '-Done'