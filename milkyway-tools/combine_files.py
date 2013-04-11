import math as ma
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import astro_coordinates as co
#import glob

""" Takes star files and puts them together, writing the total number of stars as 
the first line """

#stripes = glob.glob("./sep_lbr/Hold/*")
path = "./" #"./sep_lbr/Hold/"
prefs = ["bg-", "s2-", "s3-"]
wedges = ["09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", 
        "22", "23"] #, "79", "82", "86"]

for wedge in wedges:
    out = []
    for pref in prefs:
        try:  infile = open(path+pref+wedge+".txt", 'r')
        except IOError:  print "# - Skipping {0}{1}".format(pref, wedge);  continue
        for line in infile:
            out.append(line.strip())
        infile.close()
    length = len(out)
    if length == 0:  continue
    # Now print this to a file
    newfile = "stars-"+wedge+"-sansSgr.txt"
    outfile = open(newfile, "w")
    outfile.write(str(length))
    for item in out:
        outfile.write("\n"+item)
    outfile.close()
    print "Finished writing {0}".format(newfile)
