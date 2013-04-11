import math as ma
import numpy as np
import scipy as sc
import files as fi
import sdss_visualizers as sv

'''python script for producing plots from sdss_visualizers.py.
    Matthew Newby, November 15,2010'''

# wedges shall have:  wedge#, mu_min, mu_max, g_min, g_max
stripes = [
    #[9, 170.0, 235.0, 16.0, 23.5],
    #[10, 165.0, 245.0, 16.0, 23.0],
    #[11, 150.0, 229.0, 16.0, 23.0],
    #[12, 135.0, 235.0, 16.0, 23.0],
    #[13, 135.0, 235.0, 16.0, 23.0],
    [14, 135.0, 235.0, 16.0, 22.5],
    #[15, 135.0, 240.0, 16.0, 22.5],
    #[16, 135.0, 240.0, 16.0, 22.5],
    #[17, 135.0, 235.0, 16.0, 22.5],
    #[18, 135.0, 240.0, 16.0, 22.5],
    #[19, 135.0, 230.0, 16.0, 22.5],
    #[20, 133.0, 249.0, 16.0, 22.5],
    #[21, 133.0, 210.0, 16.0, 22.5],
    #[22, 131.0, 225.0, 16.0, 22.5],
    #[23, 133.0, 230.0, 16.0, 22.5],
    #[79, 311.0, 416.0, 16.0, 22.5],
    #[82, 310.0, 419.0, 16.0, 22.5],
    #[86, 310.0, 420.0, 16.0, 22.5]
    ]

def one_wedge_plots():
    for stripe in stripes:
        data = fi.read_data("/home/newbym2/Desktop/star_holder/stars-"+str(stripe[0])+".txt")
        wedge=stripe[0]
        mu_lim = (stripe[1], stripe[2])
        outname = "wedgeplot-"+str(stripe[0])
        sv.plot_stripe_mur(data=data, wedge=wedge, mu_lim=mu_lim, outname=outname,
                       r_lim=None, mag=0, scale=1, color=1,  nu_flatten=0)
        print "# - Stripe {0} Complete".format(stripe[0])

def multi_wedge_plots():
    prefix, pack = ["bg-", "s1-", "s2-", "s3-"], []
    folder = "/home/newbym2/Desktop/star_holder/sep_lbr/"
    for stripe in stripes:
        wedge=stripe[0]
        mu_lim = (stripe[1], stripe[2])
        outname = "multi-wedgeplot-"+str(stripe[0])
        for tag in prefix:
            pack.append(fi.read_data(folder+tag+str(wedge)+".txt"))
        sv.plot_separation_mur(pack, wedge, outname, mag=0, scale=1, color=1, mu_lim=mu_lim,
                               r_lim=None, nu_flatten=0)
        print "# - Stripe {0} Complete".format(stripe[0])

if __name__ == "__main__":
    multi_wedge_plots()

print "# - All Tasks Complete"