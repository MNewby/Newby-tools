#! /usr/bin/env python  #'Bang' line - modify as needed

import files as f
import GCplot_enc as gcp

"""This script quickly loads data and calls plotting functions.

Matthew Newby (RPI), Jan 11, 2011
"""

main_data = f.read_csv('HR_NGC_5053_cluster.csv')
back_data = f.read_csv('HR_NGC_5053_background.csv')

gcp.plot_infiles(main_data, back_data, GCname='5053', to_file=1)