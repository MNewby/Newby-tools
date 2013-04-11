import scipy as sc
import math as m
import numpy as np
import files as f
import GCplot_enc as gcp

"""Script for quick plotting of HR diagrams.

Matthew Newby, RPI, Aug 25, 2010"""

gc_name = 'NGC_6341_background'
dist = 8.2
gc_file = 'HR_' + gc_name + '.csv'
#load file
HR_data = f.read_csv(gc_file)
#plotfile
if (gcp.plot_HR_diagram(HR_data, dist, gc_name)==1):
    print '-Done'