#! /usr/bin/env python  #'Bang' line - modify as needed

import math as m
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
#import sys

"""This is a script for plotting examples of the search methods on Milkyway@home

Matthew Newby (RPI), May 02, 2013
"""

def make_surface():
    mu1, mu2 = 50.0, 50.0
    sig1, sig2 = 10.0, 10.0
    surface = sc.zeros((101,101))
    
