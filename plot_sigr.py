#! /usr/bin/env python  #'Bang' line - modify as needed

import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import sys

"""This script is for plotting the sig_r data and fit
Matthew Newby (RPI), April 13, 2010
"""

filename = "results_sigmoid_2reruns_all.txt"
data = f.read_data(filename)
runs, single = [], []
for i in range(1,len(data[:,0])):
    if a 