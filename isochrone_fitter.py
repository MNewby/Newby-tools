import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import glob

"""python script for loading and fitting An isochrones to fiducial sequences.
Matthew Newby, RPI June 2, 2011"""

#NEED:  Star file in, fiducial sequence, isochrone loading, isochrone fitting, isochrone plotting
#Metallicity, data in, vary distance, alpha as an option

class Isochrone:
    def __init__(self, age, metal, data):
        #name?
        self.age=age
        self.metal=metal
        self.data=data
    
def read_iso_file(filename):
    isofile = open(filename, "r")
    
    isofile.close()

def load_isochrones(metal, distance, alpha=0.0, type="corr"):
    files = glob.glob("../An_Isochrones/*corr*")
    #readfile = open(filename, "r")
    return 0


def two_sigma_rejection():
    return 0

def get_isochrone():  # Get (and plot?) an isocrone, given an age and metallicity
    return 0