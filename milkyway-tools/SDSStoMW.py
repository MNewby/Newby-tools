import math as ma
import numpy as np
import scipy as sc
import files as fi
import astro_coordinates as coor
import sys

def convert(filename, RaDec=0):
    readfile = open(filename, "r")
    holder, MW = [], []
    for line in readfile:
        if line[0] == "#":  continue
        if line.strip() == "": continue
        holder = line.split(",")
        for i in range(len(holder)):  holder[i] = eval(holder[i])
        if RaDec==1:
            l,b = coor.EqTolb(holder[0], holder[1])
            holder[0], holder[1] = round(l,6) , round(b,6)
        r = coor.getr(holder[2], 4.2)
        holder[2] = round(r, 6)
        MW.append(holder[:3])
    readfile.close()
    nstars = len(MW)
    out = sc.array(MW)
    fi.write_data(out, filename[:-4]+"_new.txt", " ", str(nstars))
    

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:  filename = "Stripe82_coadd.csv"
    if len(sys.argv) > 2:
        RaDec = sys.argv[2]
    else:  RaDec=0
    """ FUNCTIONS """
    convert(filename, RaDec)