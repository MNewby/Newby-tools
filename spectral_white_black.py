import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap

'''python code for generating a custom coloramp for matplotlib.
    Based on the 'spectral' color map, but with white and black reversed.
    Also contains a 
    Matthew Newby, March 19, 2012'''

cdict1 = {'red':  ((0.0, 1.0, 1.0),
                   (0.25, 0.75, 0.75),
                   (0.5, 0.5, 0.5),
                   (0.75, 0.25, 0.25),
                   (1.0, 0.0, 0.0)),
         'green': ((0.0, 1.0, 1.0),
                   (0.25, 0.75, 0.75),
                   (0.5, 0.5, 0.5),
                   (0.75, 0.25, 0.25),
                   (1.0, 0.0, 0.0)),
         'blue':  ((0.0, 1.0, 1.0),
                   (0.25, 0.75, 0.75),
                   (0.5, 0.5, 0.5),
                   (0.75, 0.25, 0.25),
                   (1.0, 0.0, 0.0))        }

white_black = LinearSegmentedColormap('grey', cdict1)

spectral_colors = fi.read_data('spectral_cm.txt')
spectral_wb = ListedColormap(spectral_colors[:,:3], name="spectral_wb", N=256)