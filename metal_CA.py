import math as ma
import numpy as np
import scipy as sc
import files as fi
import GC_CAscript_near as con
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt

'''python script for analyzing the affects of SDSS color errors on photometric
metallicity determinations.
Matthew Newby, July 21,2011'''

k=[-13.13, 14.09, 28.04, -5.51, -5.90, -58.68, 9.14, -20.61, 0.0, 58.20]


def metal_con(filename, distances, real_dist, bins=35, limits=(-3.1, 0.2),
              avgs=1, detection=1, tag="out"):
    """ main bit """
    if filename[-4:] == '.csv':  delim = ','
    else:  delim = None
    data = shift_data(fi.read_data(filename, delim), real_dist, distances[0])
    mod_actual = 5.*(ma.log10(real_dist*1000) - 1.)
    mod_new = 5.*(ma.log10(distances[0]*1000) - 1.)
    mod = mod_actual - mod_new
    print "Effective Magnitude Shift = {0}, average g={1}".format(mod, sc.mean(data[:,2]))
    new_data = cut_data(data, 4,2,3,5, deff=0, modulus=mod, full=1)
    FeH = get_photo_metal(new_data[:,4],new_data[:,2],new_data[:,3])
    ref_hist = np.histogram(FeH, bins, limits)
    hist = []
    #Also iterate over several runs and average
    for i in range(len(distances)):
        print "#- Convolving to distance {0} kpc".format(distances[i])
        if i==0:  deff=0
        else: deff=detection
        temp_hist = []
        for j in range(avgs):
            #holds dist constant, applies appropriate errors for new distance
            new_data = con.convolve(data, real_dist, distances[i])
            #shift data so detection efficiency works correctly;  has no noticable effect if deff=0
            new_data = shift_data(new_data, distances[0], distances[i])
            # apply color cuts and detection efficiency to shifted and convolved data
            new_data = cut_data(new_data, 4,2,3,5, deff=deff, modulus=None, full=0)
            print "Average g = {0}, total stars = {1}".format(sc.mean(new_data[:,2]), len(new_data[:,0]))
            FeH = get_photo_metal(new_data[:,4],new_data[:,2],new_data[:,3])
            temp_hist.append(np.histogram(FeH, bins, limits))
        new_hist = avg_hists(temp_hist)
        hist.append(new_hist)
    plot_hists(hist, ref_hist, distances, tag)
    return hist

def avg_hists(hists_in):
    """Average a set of incoming histograms """
    #If only one histogram is submitted, unpack it and return it
    if len(hists_in) == 1:
        return hists_in[0]
    new_edges = hists_in[0][1]
    # unpack the histograms from their ugly tuple
    hists = []
    for i in range(len(hists_in)):
        hists.append(hists_in[i][0])
    new_heights = []
    for i in range(len(hists[0])):  #iterate over bins
        height = 0.0
        for j in range(len(hists)):  #iterate over histograms
            height = height + hists[j][i]
        height = height / len(hists)
        new_heights.append(height)
    return (sc.array(new_heights), new_edges)

def get_photo_metal(u,g,r):
    """ Get the photometric metallicity, following appendix A1 from:
    http://adsabs.harvard.edu/abs/2010ApJ...716....1B"""
    #FeH=k[0]+k[1]*ug+k[2]*gr+k[3]*ug*gr+k[4]*ug^2+k[5]*gr^2+k[6]*ug^2*gr+k[7]*ug*gr^2+k[8]*ug^3+k[9]*gr^3
    FeH = 1000.0*sc.ones(len(u))
    for i in range(len(u)):
        ug = u[i] - g[i]
        gr = g[i] - r[i]
        part1 = k[0] + k[1]*ug + k[2]*gr + k[3]*ug*gr
        part2 = k[4]*ug*ug + k[5]*gr*gr
        part3 = k[6]*ug*ug*gr + k[7]*ug*gr*gr
        part4 = k[8]*ug*ug*ug + k[9]*gr*gr*gr
        FeH[i] = part1 + part2 + part3 + part4
    return FeH

def shift_data(data, real_dist, new_dist, ignore_cols=[0,1] ):
    """Shift data so that it appears to be at a different distance """
    mod_actual = 5.*(ma.log10(real_dist*1000) - 1.)
    mod_new = 5.*(ma.log10(new_dist*1000) - 1.)
    shift = mod_actual - mod_new
    #print "shift = {0}".format(shift)
    for i in range(len(data[:,0])):
        for j in range(len(data[0,:])):
            #skip the columns in ignore_cols
            if ignore_cols.count(j) > 0:  continue
            # otherwise, shift data
            else:  data[i,j] = (data[i,j] - shift)
    return data
    
def clip_data():
    """ Remove data outside a certain distance range, only need this for non-clusters """
"""
0.2 < g-r < 0.6
0.7 < u-g < 2.0
-0.25 < (g-r)-0.5(u-g) < 0.05
-0.2 < 0.35(g-r)-(r-i) < 0.10
"""

def cut_data(data, u_col, g_col, r_col, i_col, deff=0, modulus=None, full=0):
    """ Remove data outside a certain color range """
    ug_range = [0.7, 2.0]
    gr_range = [0.2, 0.6]
    grug_range = [-0.25, 0.05]
    grri_range = [-0.2, 0.10]    
    new_data = []
    kills = 0
    for i in range(len(data[:,0])):
        ug = data[i,u_col] - data[i,g_col]
        gr = data[i,g_col] - data[i,r_col]
        if (ug < ug_range[0]) or (ug > ug_range[1]):  continue
        if (gr < gr_range[0]) or (gr > gr_range[1]):  continue
        if full == 1:
            ri = data[i,r_col] - data[i,i_col]
            grug = gr - 0.5*ug
            if (grug < grug_range[0]) or (grug > grug_range[1]):  continue
            grri = 0.35*gr - ri
            if (grri < grri_range[0]) or (grri > grri_range[1]):  continue
        if (deff==1):
            if (np.random.uniform() > sigmoid_error(data[i,g_col], modulus) ):
                kills = kills + 1
                #print data[i,g_col]
                continue
        new_data.append(data[i,:])
    print "keepers: {0}, color kills: {1}, deff kills: {2}".format(
                        len(new_data), (len(data[:,0]) - len(new_data)- kills), kills)
    return sc.array(new_data)

def sigmoid_error(x, modulus=None):
    """Application of detection efficiency"""    
    s = [0.9402, 1.6171, 23.5877]
    if modulus != None:
        s[2] = s[2] - modulus
    detection_efficiency = s[0] / (np.exp(s[1]*(x - s[2])) + 1.)
    return detection_efficiency

def plot_hists(histograms, ref, ids=None, tag="out"):
    """ Ideally, plot the ref with each other hist in a subplot.  maybe also plot
        the detection-efficiency-subtracted data with other data"""
    o_stars = sc.sum(ref[0])
    for i in range(len(histograms)):
        t_stars = sc.sum(histograms[i][0])
        fig = plt.figure(i+1)
        plt.bar(histograms[i][1][:-1], histograms[i][0], 0.1)
        plt.bar(ref[1][:-1], ref[0], 0.1, fill=False, edgecolor='r', ls="dashed")
        plt.xlabel(r"[Fe/H]")
        plt.ylabel("counts")
        if ids != None:
            out = "{0} kpc, {1} stars,\n{2} of original".format(ids[i],
                            t_stars, round( (float(t_stars)/float(o_stars)), 2 ) )
            plt.text(0.0, 120.0, out, fontsize=12)
        plt.savefig(tag+"_"+str(ids[i])[:5]+".ps", papertype='letter')
    #plt.show()

if __name__ == "__main__":
    d = sc.array([7.7, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0])
    #d = sc.arange(1.0, 82.0, 10.0)
    tag = "real_noFx"
    averages = 10
    limits = (-4.0, 2.0)
    bins = (limits[1] - limits[0]) / 0.1
    print "limits = {0}, bins = {1}, averages = {2}".format(limits, bins, averages)
    output = metal_con("NGC_6205_nofield.txt", d, 7.7, bins=bins, limits=limits,
              avgs=averages, detection=1, tag=tag)
    for i in range(len(output)):
        out = sc.zeros( (len(output[i][0]),2), float)
        out[:,0], out[:,1] = output[i][0], output[i][1][:-1]
        #print out.shape
        fi.write_data(out, "hist_"+tag+"_"+str(d[i])[:5]+".txt", header=tag+"_"+str(d[i])[:5])
    #do background as well?
    """ Address distance, color cuts? Distances only make sense when convolving a
    cluster, need to modify for convolution of background."""
    print "# --- Done"