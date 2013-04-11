import math as ma
import numpy as np
import scipy as sc
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import astro_coordinates as co


def load_params(filename):
    """ loads a set of best-fit parameters from a file """
    holder = []
    readfile = open(filename, "r")
    for line in readfile:
        if line.strip() == "":  continue
        if line[:2] == "ps":  continue
        if line[:2] == "de":  continue
        temp = line.split(',')[2:]    
        for i in range(len(temp)):
            temp[i] = temp[i].strip("[] \n")
        #print temp, "\n"
        holder.append(temp)
    return holder
    

stripes = ['09']+range(10,23)
paramstr = load_params("../sgrOther/Results_de.txt")
print "# <stripe> -s<stream#>:  <weight>; <l>, <b>, <r>"
for i in range(len(stripes)):
    #print paramstr[i]
    l,b,r = co.GC2lbr(float(paramstr[i][3]), 0.0, float(paramstr[i][4]), int(stripes[i]))
    print "{0} -s1:  {1}; {2}, {3}, {4}".format(int(stripes[i]),paramstr[i][2],l,b,r)
    l,b,r = co.GC2lbr(float(paramstr[i][9]), 0.0, float(paramstr[i][10]), int(stripes[i]))
    print "{0} -s2:  {1}; {2}, {3}, {4}".format(int(stripes[i]),paramstr[i][8],l,b,r)
    l,b,r = co.GC2lbr(float(paramstr[i][15]), 0.0, float(paramstr[i][16]), int(stripes[i]))
    print "{0} -s3:  {1}; {2}, {3}, {4}\n".format(int(stripes[i]),paramstr[i][4],l,b,r)
    #print "\n"
print "# --- Done"
