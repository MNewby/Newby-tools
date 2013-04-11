import math as ma
import numpy as np
import scipy as sc
import files as fi
import astro_coordinates as coor
import test_wedge_generator as twg

# ditch params class?

# keep stars that break mu, r limits and detection efficiency - record separately;
# kill all stars that break nu limits, except in southern stripes?

# generate stream consistent with Separation results

def integrate_stripes(stripes, fit):
    mu_min = [170.0, 165.0, 150.0, 135.0, 135.0, 135.0, 135.0, 135.0, 135.0, 135.0,
              135.0, 133.0, 133.0, 131.0, 133.0, 311.0, 310.0, 310.0]
    mu_max = [235.0, 227.0, 229.0, 235.0, 235.0, 235.0, 240.0, 240.0, 235.0, 240.0,
              230.0, 249.0, 210.0, 225.0, 230.0, 416.0, 419.0, 420.0]
    N_stars = [19386, 18734, 11096, 17044, 18736, 15409, 12519, 12248, 8853, 7328,
               5479, 4450, 3486, 2425, 971, 9511, 16119, 16603]
    deg = 180.0/ma.pi
    for i in range(len(stripes)):
        if stripes[i] != 82:  continue
        params = twg.ParamSet()
        params.wedge = stripes[i]
        params.background = [0.0, 1.0, float(fit[i][1]), float(fit[i][2]), 1.0]
        # flip the weird angles, and convert to degrees
        if stripes[i]==9 or stripes[i]==22 or stripes[i]==82:
            theta, phi = (ma.pi-eval(fit[i][6]) ), (eval(fit[i][7]) + ma.pi)
        else:
            theta, phi = eval(fit[i][6]), eval(fit[i][7])
        theta, phi = coor.angle_bounds3(theta*deg, phi*deg, phi_max=180.0)
        params.streams = [ [float(fit[i][3]), float(fit[i][4]), float(fit[i][5]),
                            94.682, 116.0, float(fit[i][8])] ]  #normal to 82
        #params.streams = [ [float(fit[i][3]), float(fit[i][4]), float(fit[i][5]),
        #                    theta, phi, float(fit[i][8])] ]
        # Set g limits
        if (stripes[i]==9) or (stripes[i]==10):  g_max = 23.5
        elif (stripes[i]==11) or (stripes[i]==12):  g_max = 23.0
        else:  g_max = 22.5
        params.stripe = [(mu_min[i], mu_max[i], 10),(-1.25, 1.25, 10),(16.0,g_max,10)]
        params.update_refs()
        params.print_params()
        name = "sim_stripe_"+str(stripes[i])+"_NormNoCon.txt"
        print "STARTING: "+name
        stream_gen(params, 0, N_stars[i], fileout=name, convolve=0)

def stream_gen(params, sn, N_stars, batch=1000, fileout="streamgen82.txt",
                       detection=1, convolve=1, append=0):
    """ sn is stream number"""
    # Initialize file
    if append==0:
        out = open(fileout, 'w')
        out.write("# "+str(N_stars)+" stars, l,b,r,flag \n")
        out.close()
    # Get constants
    mu,R,theta,phi,sigma,wedge = \
        params.mu[sn],params.R[sn],params.theta[sn],params.phi[sn],params.sigma[sn], params.wedge
    u_min, u_max = twg.get_stream_length(params, sn, accuracy=0.00001)
    nu_min, nu_max = params.nu_lim[0], params.nu_lim[1]
    mu_min, mu_max = params.mu_lim[0], params.mu_lim[1]
    g_min, g_max = params.stripe[2][0], params.stripe[2][1]
    print "# - Generating Stream {0}, using parameters {1}, {2}, {3}, {4}, {5}".format(
        sn, mu, R, theta, phi, sigma)
    X_out, Y_out = 0, 0 # X_out is stars detected in stripe, Y_out is stars generated that should be there.
    while X_out < N_stars:
        mu_killed, nu_killed, mu_saved = 0,0,0
        u,v,w = twg.generate_stream(batch, u_min, u_max, sigma)
        holder = []
        for i in range(len(u)):
            mu1, nu1, r1 = coor.streamToGC(u[i],v[i],w[i],mu,R,theta,phi,wedge)
            l,b,r1 = coor.GC2lbr(mu1, nu1, r1, wedge)
            # Convolve
            if convolve==1:
                r1 = twg.star_convolution(r1)
            # Test nu, kill if invalid
            if (nu1 < nu_min) or (nu1 > nu_max):
                nu_killed=nu_killed+1
                #Y_out=Y_out+1   # COMMENT
                #holder.append([round(l,6),round(b,6),round(r1,6), 1])  # COMMENT
                continue
            # Test mu, keep in Y if out of bounds
            if (mu_max > 360.0):
                if (mu1 > (mu_max-360.0)) and (mu1 < mu_min):
                    mu_killed=mu_killed+1
                    Y_out=Y_out+1
                    holder.append([round(l,6),round(b,6),round(r1,6), 1])
                    continue
            else:
                if (mu1 < mu_min) or (mu1 > mu_max):
                    mu_killed=mu_killed+1
                    Y_out=Y_out+1
                    holder.append([round(l,6),round(b,6),round(r1,6), 1])
                    continue
            # Detection
            if detection==1:
                m_g = coor.getg(r1)
                if np.random.uniform() > (twg.sigmoid_error(m_g)):
                    Y_out=Y_out+1
                    holder.append([round(l,6),round(b,6),round(r1,6), 1])
                    continue
            # test g-limits
            if (coor.getg(r1) < g_min) or (coor.getg(r1) > g_max):
                Y_out=Y_out+1
                holder.append([round(l,6),round(b,6),round(r1,6), 1])
                continue
            # When a point passes all the tests add it to the set
            #co.stream2xyz(u[i],v[i],w[i],mu,R,theta,phi,wedge)
            holder.append([round(l,6),round(b,6),round(r1,6), 0])
            X_out = X_out + 1
        if X_out > N_stars:
            slice = -1*(X_out-N_stars)
            holder = holder[:slice]
            print "#---Sliced {0} stars to make quota".format(str(-1*slice))
            X_out = X_out + slice #Slice is negative
        #append X and Y to files
        if len(holder) != 0:
            if fi.append_data(sc.array(holder), fileout, delimiter=" ") == 1:
                print "#---Stream Progress: {0} stars of stars generated out of {1} detected".format(Y_out+X_out, N_stars)
                print "# !!! - out of mu: {0}, nu_killed: {1}".format(mu_killed, nu_killed)
    print "#---Stream {0} generation succeeded, written as {1}".format(sn, fileout)
    print " ### -----  Stream {0} Complete ----- ### \n"
    return fileout

