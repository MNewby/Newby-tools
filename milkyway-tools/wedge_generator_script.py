import sys
sys.path.insert(0, '../utilities')
import math as ma
import numpy as np
import scipy as sc
import files as fi
import test_wedge_generator as twg

""" Script for easy access to the defs in "test_wedge_generator.py"

Matthew Newby, Feb 27, 2012 """

if __name__ == "__main__":
    """
    params = twg.ParamSet()
    params.wedge = 14
    params.background = [0.0, 1.0, 0.56, 8.0, 1.0]
    params.streams = [ [-1.0, 183.946, 28.507, 1.853, -0.213, 4.0],
                       [-2.0, 213.567, 24.957, 1.834, -2.255, 2.0]
                       ]
    params.stripe = [(135.0, 235.0, 10), (-1.25, 1.25, 10), (16.0, 22.5, 10)]
    params.update_refs()
    twg.build_stripe(params, filename="stars-t14-20-3.txt", num_stars=110000,
        perturb_weight=0.20, perturb_params=(0.0, 0.79, -19.9), con=1, det=1, app=1)
    params.print_params()
    print "\n ### Done"
    """
    # Make naked streams
    stripes = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    nStars = [19386, 18734, 11095, 17044, 18736, 15409, 12519, 12248, 8853, 7328, 5479, 4450, 3486, 2425, 971]
    #read params: "/home/newbym2/Dropbox/Research/sgrnorth_paper/results_MW_easyread.txt"
    paramfile = open("/home/newbym2/Dropbox/Research/sgrnorth_paper/results_BG_easyread.txt", "r")
    p = []
    for line in paramfile:
        if line.strip() == "":  continue
        if line.strip()[0] == "#":  continue
        if line.strip()[:3] == "---":  continue
        p.append(line.split(','))
    paramfile.close()
    # loop over stripes
    # build streams
    for i in range(len(stripes)):
        params = twg.ParamSet()
        params.wedge = stripes[i]
        params.streams = [ p[i][3:9] ]
        print params.streams
        for j in range(len(params.streams[0])):
            #if j==1:  params.streams[0][j] = float(params.streams[0][j]) + 15.0
            #elif j==5:  params.streams[0][j] = float(params.streams[0][j])*2.0
            #else:     params.streams[0][j] = float(params.streams[0][j])
            params.streams[0][j] = float(params.streams[0][j])
        if stripes[i] > 12:  outerg = 22.5
        elif stripes[i] > 10:  outerg = 23.0
        else:  outerg = 23.5
        params.stripe = [(135.0, 240.0, 10), (-1.25, 1.25, 10), (16.0, outerg, 10)]
        params.update_refs()
        twg.stream_into_stripe(params, 0, int(nStars[i]), batch=1000, 
            fileout="streamgen_sgr_sim.txt", detection=1, convolve=1, append=1)
    
#if MAKE theta, phi BE IN DEGREES!!!!  <- What the hell does that mean???
