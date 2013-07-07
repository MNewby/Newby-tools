import sys
sys.path.insert(0, '../utilities')
import math as ma
import numpy as np
import scipy as sc
import files as fi
import sdss_visualizers as sdss
import sys

    
if __name__ == "__main__":
    if len(sys.argv) > 1: data = fi.read_data(sys.argv[1])
    if len(sys.argv) > 2: wedge = int(sys.argv[2])
    else:  wedge = 82
    #plot_stripe_lb(data)
    #plot_stripe_3D(data)
    #sdss.plot_stripe_mur(data, wedge, mag=0, scale=1, r_lim=(0.0, 46.0), 
    #    vm=10.0, mu_lim=(310.0, 419.0), color=1)
    #wedge = 10
    #data0 = fi.read_data("../../sep_lbr/bg-10.txt")
    #data1 = fi.read_data("../../sep_lbr/s1-10.txt")
    #data2 = fi.read_data("../../sep_lbr/s2-10.txt")
    #data3 = fi.read_data("../../sep_lbr/s3-10.txt")
    #sdss.plot_separation_mur(wedge, data0, data1, data2, data3,
    #                    outname=None, mag=0, scale=1, color=0, mu_lim=(165.0, 225.0), 
    #                    r_lim=(0.0, 70.0), vm=10.0, nu_flatten=0, bar=0)
    sdss.plot_stripe_mur(data, wedge, mag=0, scale=1, r_lim=(0.0, 55.0), vm=10.0, mu_lim=(165.0, 245.0))
    #plot_stripe_mug(data, wedge)
    #plot_stripe_mu(data, wedge)
    #plot_stripe_munu(data, wedge)
    #plot_wedge_density(data, wedge, q=0.458, r0=19.4, name=sys.argv[1][:-4]+"_rho", mag=0)
    #plot_double_mur(data, wedge=82, mu_min=310.0, mu_max=419.0,
    #               name=sys.argv[1][:-4]+"_mur", mag=0, rsize=1.0, musize=1.0)
    #do_wedgeplots()
    #do_testwedgeplots()
    #do_densityplots()
    #do_onewedgeplot("stripe-10-DR7-clean.txt", 10, None)
    #do_compare_wedges(file1="unmatched-82-1arcsec.txt", file2="matched-82-1arcsec.txt", stripe=82, hern=0)
    sdss.lbpolar_plot("/home/newbym2/Dropbox/Research/sep_lbr_sSgr/", hemi='N', bin_size=0.5, 
        outfile=None, infile="/home/newbym2/Dropbox/Research/sgrOther/sanSgr_hist0.txt", color=1, scale=1, primary=1)
    #sdss.lbpolar_plot("/home/newbym2/sscon_lbr/", hemi='N', bin_size=0.5, 
    #    outfile=None, infile=None, color=1, scale=1, primary=1)
    #sdss.sgr_xz_plot(folder="./sep_lbr/", bin_size=0.5, outfile=None, infile=None,
    #            primary=1, color=1, scale=0)
    # /home/newbym2/Desktop/star_holder/sep_lbr
    #sdss.sgr_plot()
    """t0 = time.time()
    wedges = range(70,90)
    all_stripes(wedges, "FTO_south_mnewby.csv")
    tf = time.time()
    print "Total time elapsed: {0} minutes".format(((tf-t0)/60.0))"""
    
    print "# --- Done"
