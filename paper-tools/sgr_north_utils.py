import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import astro_coordinates as coor
import sgr_law as sl
import make_latex_table as tex
import sgr_density_analysis as sda

rad = ma.pi / 180.0
deg = 180.0 / ma.pi
sun = [-8.5, 0.0, 0.0]
GC = [0.0, 0.0, 0.0]
Sgr_lbr = [5.6, -14.2, 24.0]
Sgr = coor.lbr2xyz(Sgr_lbr[0], Sgr_lbr[1], Sgr_lbr[2])

'''python scripts for producing components of the sgrNorth paper.
Matthew Newby, April 16, 2012'''

def load_params(filename):
    stripes, params = [], []
    readfile = open(filename, "r")
    for line in readfile:
        if line.strip() == "":  continue
        if line.strip()[0] =="#":  continue
        if line.strip()[0:3] == "---":
            stripes.append(int(line.strip().strip(" -Stripe")))
            #print "Loaded stripe {0}".format(stripes[-1])
            continue
        params.append(line.strip().split(","))
    readfile.close()
    return stripes, params

def results_to_table(stripes, params):
    """ Build latex tables from a file of results """
    # Background params
    bg_data = []
    for i, stripe in enumerate(stripes):
        #print i, stripe, params[i]
        bg_data.append([stripe, params[i][1], params[i][2] ])
    bg_table = tex.DataTable(sc.array(bg_data), 'lcc', "Spheroid Parameters", "spheroidtable",
                          ["Stripe","$q$","$r_0$"])
    tex.deluxe_table(bg_table, roundoff=3)
    # Sgr params
    sgr_data = []
    for i, stripe in enumerate(stripes):
        if stripe==9 or stripe==22 or stripe==82:
            theta, phi = (ma.pi-eval(params[i][6]) ), (eval(params[i][7]) + ma.pi) #not sure if needed anymore
        else:
            theta, phi = eval(params[i][6]), eval(params[i][7])
        theta, phi = coor.angle_bounds3(theta*deg, phi*deg, phi_max=180.0)
        theta, phi = theta*rad, phi*rad
        sgr_data.append([stripe, params[i][3], params[i][4], params[i][5], str(theta),
                         str(phi), params[i][8] ])
    sgr_table = tex.DataTable(sc.array(sgr_data), 'lcccccc', "Sagittarius Stream Parameters",
                              "sgrtable", ["Stripe","$\epsilon$","$\mu$ ({}^{\circ})", "R (kpc)",
                                           "\theta (rad)", "\phi (rad)", "\sigma"])
    tex.deluxe_table(sgr_table, roundoff=3)

def make_table_werrors(stripes, params):
    file2="../sgrnorth_paper/results_errors_final.txt"
    errors=fi.read_data(file2)
    for i, stripe in enumerate(stripes):
        if stripe==9 or stripe==22 or stripe==82:
            theta, phi = (ma.pi-eval(params[i][6]) ), (eval(params[i][7]) + ma.pi) #not sure if needed anymore
        else:
            theta, phi = float(params[i][6]), float(params[i][7])
        theta, phi = coor.angle_bounds3(theta*deg, phi*deg, phi_max=180.0)
        theta, phi = theta*rad, phi*rad
        fit = [stripe, float(params[i][3]), float(params[i][4]), float(params[i][5]),
               theta, phi, float(params[i][8]) ]
        # put it together
        out = str(stripe)+" & "
        for j in range(2,8):
            out = out + str(round(fit[j-1],2))+r" \pm "+str(round(errors[i,j],2))+r" & "
        out = out + r" \\"
        print out, "\n"
        back = [stripe, ]

def plane_fit_table():
    """  Make latex tables from a set of plane-fits """
    data = fi.read_data("../sgrnorth_paper/planefit_results.txt", ",")
    for i in range(len(data[:,0])):
        string = ""
        for j in range(4):
            string = string + " & " + str(round(data[i,j],3)) + r"$\pm$" + str(round(data[i,j+4],3))
        print string

def sgr_xyz_plots(stripes, params):
    """ Plots stream positions & directions in galactic XYZ"""
    positions = sc.zeros((len(stripes),3))
    vectors = sc.zeros((len(stripes),3))
    for i, stripe in enumerate(stripes):
        positions[i,0], positions[i,1], positions[i,2] = coor.GC2xyz(
            eval(params[i][4]), 0.0, eval(params[i][5]), stripe)
        if stripe==9 or stripe==22 or stripe==82:
            theta, phi = (ma.pi-eval(params[i][6]) ), (eval(params[i][7]) + ma.pi)
        else:
            theta, phi = eval(params[i][6]), eval(params[i][7])
        theta, phi = coor.angle_bounds3(theta*deg, phi*deg, phi_max=180.0)
        theta, phi = theta*rad, phi*rad
        print stripe, theta, phi
        vectors[i,0] = ma.sin(theta)*ma.cos(phi)
        vectors[i,1] = ma.sin(theta)*ma.sin(phi)
        vectors[i,2] = ma.cos(theta)
    off = 1.0
    """ Plot xy """
    plt.figure(1)
    #ax = plt.subplot(111)
    plt.scatter(GC[0], GC[1], s=40, edgecolor="k", marker="o", facecolor="black")
    plt.text(GC[0]+off, GC[1]+off, "GC", fontsize=10)
    plt.scatter(sun[0], sun[1], s=20, edgecolor="k", marker=(8,2,0), facecolor="white")
    plt.text(sun[0]+off, sun[1]+off, "Sun", fontsize=10)
    plt.scatter(Sgr[0], Sgr[1], s=20, edgecolor="k", marker="o", facecolor="black")
    plt.text(Sgr[0]+off, Sgr[1]+off, "Sgr Core", fontsize=10)
    t = sc.arange(0.0, 2.0*ma.pi, 0.01)
    plt.plot(30.0*sc.cos(t), 30.0*sc.sin(t), "k:")
    offset_x = [0.0, 0.0, 0.0, -1.0, -2.0, -1.0, -1.0, 0.0, 0.0, 0.0, -1.0,
                -2.0, -3.2, -3.0, 0.5, 0.0, 0.0, -1.0]
    offset_y = [0.0, 0.0, -1.0, -1.9, 0.2, 0.5, 0.5, -1.0, -2.0, -3.0, -3.0,
              -2.9, -2.5, 0.5, 0.5, -2.0, -2.0, -2.25]
    for i in range(len(stripes)):
        #plt.annotate(str(stripes[i]), xy=(positions[i,0], positions[i,1]),
        #            xytext=( (positions[i,0]+(vectors[i,0]*100.0)), (positions[i,1]+(vectors[i,1]*100.0)) ),
        #            arrowprops=dict(facecolor='black', arrowstyle="->") )
        if stripes[i] < 75:  color='white'
        else:  color='black'
        plt.arrow(positions[i,0], positions[i,1], vectors[i,0]*10.0, vectors[i,1]*10.0,
                  width=0.20, head_width=1.5, facecolor=color)
        if stripes[i] < 20 or stripes[i] > 19:
            plt.text(positions[i,0]+offset_x[i], positions[i,1]+offset_y[i],
                     r"$"+str(stripes[i])+r"$", fontsize=10)
    plt.xlabel(r"$X_{GC}$ (kpc)")
    plt.ylabel(r"$Y_{GC}$ (kpc)")
    plt.axis('equal')
    plt.ylim(-30, 30)
    plt.xlim(-30, 30)
    locs, labels = plt.xticks()
    for i in range(len(labels)):  labels[i] = r"$"+str(locs[i])+r"$"
    plt.xticks(locs, labels, fontsize=12)
    locs, labels = plt.yticks()
    for i in range(len(labels)):  labels[i] = r"$"+str(locs[i])+r"$"
    plt.yticks(locs, labels, fontsize=12)
    """ Plot xz """
    plt.figure(2)
    #ax = plt.subplot(111)
    plt.scatter(GC[0], GC[2], s=40, edgecolor="k", marker="o", facecolor="black")
    plt.text(GC[0]+off, GC[2]+off, r"GC", fontsize=10)
    plt.scatter(sun[0], sun[2], s=20, edgecolor="k", marker=(8,2,0), facecolor="white")
    plt.text(sun[0]+off, sun[2]+off, "Sun", fontsize=10)
    plt.scatter(Sgr[0], Sgr[2], s=20, edgecolor="k", marker="o", facecolor="black")
    plt.text(Sgr[0]+off, Sgr[2]+off, "Sgr Core", fontsize=10)
    plt.plot([-30.0, 30.0], [0.0, 0.0], "k:")
    offset_x = [0.2, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, -0.5, 1.0, 1.0, 0.5,
                0.75, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0]
    offset_z = [0.0, 0.0, 0.0, -2.5, 0.0, 0.0, -2.0, 0.5, -1.0, -1.0, 0.25,
                -1.0, -0.5, 0.0, -0.5, 0.0, 0.5, 0.5]
    for i in range(len(stripes)):
        #plt.annotate(str(stripes[i]), xytext=(positions[i,0], positions[i,2]),
        #            xy=( (positions[i,0]+(vectors[i,0]*10.0)), (positions[i,2]+(vectors[i,2]*10.0)) ),
        #            arrowprops=dict(facecolor='black', arrowstyle="->") )
        if stripes[i] < 75:  color='white'
        else:  color='black'
        plt.arrow(positions[i,0], positions[i,2], vectors[i,0]*10.0, vectors[i,2]*10.0,
                  width=0.25, head_width=2.0, facecolor=color)
        if stripes[i] < 20 or stripes[i] > 19:
            plt.text(positions[i,0]+offset_x[i], positions[i,2]+offset_z[i],
                     r"$"+str(stripes[i])+r"$", fontsize=10)
    plt.xlabel(r"$X_{GC}$ (kpc)")
    plt.ylabel(r"$Z_{GC}$ (kpc)")
    plt.axis('equal')
    plt.xlim(-45, 45)
    plt.ylim(-40, 50)
    locs, labels = plt.xticks()
    for i in range(len(labels)):  labels[i] = r"$"+str(locs[i])+r"$"
    plt.xticks(locs, labels, fontsize=12)
    locs, labels = sc.arange(-40.0, 50.0, 20.0), []
    for i in range(len(locs)):  labels.append(r"$"+str(locs[i])+r"$")
    plt.yticks(locs, labels, fontsize=12)
    # done
    plt.show()

def get_errors(file1="../sgrnorth_paper/results_BG_easyread.txt",
               file2="../sgrnorth_paper/results_errors.txt"):
    """ Takes in mu,nu,r errors, and turns them into a 3D "sigma sphere" for
        fitting a plane - GARBAGE"""
    wedges = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 79, 82, 86]
    data = fi.read_data(file1, ",")
    errors = fi.read_data(file2)
    circ = []  # circular errors
    for i in range(len(data[:,0])):
        # Get x,y,z errors due to d_R
        x0, y0, z0 = coor.stream2xyz(0.0, 0.0, 0.0, data[i,4], errors[i,4], 0.0, 0.0, wedges[i])
        # Get x,y,z errors due to d_mu - the difference between the original vector, and that vector shifted by d_mu
        x1, y1, z1 = coor.stream2xyz(0.0, 0.0, 0.0, data[i,4], data[i,5], 0.0, 0.0, wedges[i])
        x2, y2, z2 = coor.stream2xyz(0.0, 0.0, 0.0, (data[i,4]+errors[i,3]), data[i,5], 0.0, 0.0, wedges[i])
        x3, y3, z3 = (x2-x1), (y2-y1), (z2-z1)
        #print x0,y0,z0, x3,y3,z3
        X = sc.sqrt(x0*x0 + x3+x3)
        Y = sc.sqrt(y0*y0 + y3+y3)
        Z = sc.sqrt(z0*z0 + z3+z3)
        print X,Y,Z
        circ.append(float(sc.sqrt(X*X+Y*Y+Z*Z) ) )
    for item in circ:  print item
    return circ

def get_distances():
    north = [-0.203, 0.933, 0.298, -1.193]
    south = [0.0326, 0.988, 0.153, -1.201]
    print plane_distance(sc.array(Sgr), north)
    print plane_distance(sc.array(Sgr), south)
    print plane_distance(sc.array([0.0,0.0,0.0]), north)
    print plane_distance(sc.array([0.0,0.0,0.0]), south)    

def plane_distance(data, params):
    a, b, c, d = params
    x, y, z = data[0], data[1], data[2]
    # Plane normal:  V = [a,b,c]
    bottom = np.sqrt(a*a + b*b + c*c)
    # Distance, Wolfram Mathworld
    top = a*x + b*y + c*z + d
    return top / bottom

def make_coord_table(stripes, params):
    """ Makes a table of the Sgr dections, in standard coordinates
        Stripe #, x, y, z, l, b, lam, beta, new_lam, new_beta, R"""
    Sgr_lbr = [5.6, -14.2, 24.0]
    Sgr_xyz = [14.655645800014774, 2.2704309216189817, -5.8873772610912614]
    Stream_North = [-0.199, 0.935, 0.293, -0.988]
    Stream_South = [0.032, 0.987, 0.158, -1.073]
    coords = []
    for i, stripe in enumerate(stripes):
        mu, r, theta, phi = eval(params[i][4]), eval(params[i][5]), \
            eval(params[i][6]), eval(params[i][7])
        x, y, z = coor.stream2xyz(0.0, 0.0, 0.0, mu, r, theta*deg, phi*deg, stripe)
        l,b,r = coor.xyz2lbr(x,y,z)
        Xs, Ys, Zs, lam, beta, rs = sl.law_xyz2sgr_sun(x,y,z)
        if lam > 180.0:  lam = lam - 360.0
        # Just going to use Sgr North plane for now
        xp, yp, zp = coor.xyz2plane(x,y,z, new_x=Sgr_xyz, plane=Stream_North, origin=-8.5)
        long, lat, rp = coor.xyz2longlat(xp,yp,zp)
        print r, rs, rp
        coords.append([stripe, x, y, z, l, b, lam, beta, long, lat, r])
    # Now make table
    coord_table = tex.DataTable(sc.array(coords), 'lcccccccccc', "Sagittarius Stream Centers",
            "coordtable", ["Stripe",r"$X_{GC}$",r"$Y_{GC}$",r"$Z_{GC}$",r"$l$",r"$b$",
                           r"$\Lambda_{LM}$",r"$\Beta_{LM}$",r"$\Lambda_{new}$",
                           r"$\Beta_{new}$", r"$R_{\sun} (kpc)$"])
    tex.deluxe_table(coord_table, roundoff=1)
    
    """Sgr_lbr = [5.6, -14.2, 24.0]
    Sgr = [14.655645800014774, 2.2704309216189817, -5.8873772610912614]
    Law = [-0.064, 0.970, 0.233, 0.23]  #GC-centered
    Law_sun = [-0.064, 0.970, 0.233, 0.776]  #Sun-centered
    New_North = [-0.203, 0.933, 0.298, -1.193]
    l0 = [229.1, 217.2, 209.0, 204.0, 201.8]
    b0 = [75.6, 72.9, 67.5, 62.4, 53.2]
    g0 = [21.75, 21.45, 21.25, 21.65, 21.40]
    #lam = [75.0, 85.0, 95.0, 105.0, 115.0]  
    #beta = []  # Secondary peak, Koposov 2012
    beta = [-1.3, -2.5, -1.3, -0.9, -1.4]  # Primary peak, Koposov 2012
    #print "POLE, Law:", coor.xyz2lbr(-Law[0], Law[1], Law[2], d0=0.0)
    #print "POLE, New_North:", coor.xyz2lbr(New_North[0], New_North[1], New_North[2], d0=0.0)
    for i in range(len(lam)):
        x1,y1,z1 = coor.lbr2xyz(lam[i], beta[i], 1.0, d0=0.0)
        print x1,y1,z1
        x,y,z = coor.plane2xyz(x1,y1,z1, Sgr, Law)
        print x,y,z
        l,b,r = coor.xyz2lbr(x,y,z)
        print l,b,r
        ra, dec = coor.lbToEq(l,b)
        print ra, dec
        print
    #r = coor.getr(21.45)
    #x,y,z = coor.lbr2xyz(217.2, 72.9, r)
    #lam, beta, r = coor.xyz2lambeta(x+8.5,y,z, Sgr, Law_sun)
    #print lam, beta, r
    #x, y, z = coor.lambeta2xyz(l,b,r, Sgr, Law)
    #print x, y, z
    #print coor.xyz2lambeta(0.0, 3.0, 1.0, [15.0, 0.0, 0.0], [0.0, 1.0, 0.0, -2.0])
    #print sl.law_xyz2sgr_sun(x,y,z)
    #print sl.law_xyz2sgr_gal(10.0, 0.0, 0.0)
    #x,y,z = coor.rot_3D(Sgr[0], Sgr[1], Sgr[2], "zxz", [-3.77*rad, -76.5*rad, -22.2*rad])
    #print coor.xyz2lbr(x,y,z,d0=0.0)"""
    
def astro_test():
    """ """
    Sgr_lbr = [5.6, -14.2, 24.0]
    Sgr = coor.lbr2xyz(Sgr_lbr[0], Sgr_lbr[1], Sgr_lbr[2])
    print Sgr
    # (14.655645800014774, 2.2704309216189817, -5.8873772610912614)
    Law = [-0.064, 0.970, 0.233, 0.23]  #GC-centered
    Mine = [-0.203, 0.933, 0.298, -1.193]
    lam = [75.0, 85.0, 95.0, 105.0, 115.0]
    beta = [9.5, 10.3, 10.3, 10.7, 11.8]
    
    t = sc.arange(0.0, 2.0*ma.pi, ma.pi/16.0)
    xx = 10.0*sc.cos(t)
    zz = 10.0*sc.sin(t)
    
    x = xx #[Sgr[0], 1.0, 0.0]
    y = sc.zeros(len(xx)) #[Sgr[1], 0.0, 1.0]
    z = zz #[Sgr[2], 0.0, 0.0]
    # XXXXXXXX  NEED SUN-CENTER LAW PLANE!
    lps, bps, rps = [], [], []
    lams, betas, rs = [], [], []
    for i in range(len(x)):
        xp, yp, zp = coor.xyz2plane(-x[i],y[i],z[i], new_x=Sgr, plane=Law, origin=sun[0])
        lp, bp, rp = coor.xyz2longlat(xp,yp,zp)
        lps.append(lp); bps.append(bp); rps.append(rp)
        print lp, bp, rp
        Xs, Ys, Zs, lam, beta, r = sl.law_xyz2sgr_sun(x[i],y[i],z[i]) #, Xsun=dsun)
        print lam, beta, r
        lams.append(lam); betas.append(beta); rs.append(r)
    plt.figure(1)
    plt.scatter(lps, bps, c='b')
    plt.scatter(lams, betas, c='g')
    plt.show()
    
    #print Sgr
    #print coor.xyz2plane(2.0, 0.0, 0.0, [0.0,1.0,0.0], [0.0, 0.0, 1.0, 0.0])
    #print xyz2plane(Sgr[0], Sgr[1], Sgr[2], [2.22222222, 3.333333333, 4.444444444], [0.0326, 0.988, 0.153, 0.0])
    
    #print coor.xyz2plane(Sgr[0], Sgr[1], Sgr[2], list(Sgr), Law)
    #print coor.xyz2lambeta(Sgr[0], Sgr[1], Sgr[2], list(Sgr), Law)
    #print coor.xyz2plane_old(Sgr[0], Sgr[1], Sgr[2], list(Sgr), Law)
    #print coor.plane2xyz(15.955954308821331, 4.1893571944839891e-16, 0.12264161315188961, list(Sgr), Law)
    
    #print plane2xyz(24.0-8.5, 0.0, 0.0, list(Sgr), Law)
    #print plane2xyz(15.955659255114893, -0.068783217539582636, 0.14044740696010716, list(Sgr), Law)
    #print coor.rot_3D(Sgr[0], Sgr[1], Sgr[2], "zxz", [-183.8*rad, 76.5*rad, 201.6*rad])
    print
    
def dotties():
    data = fi.read_data("../sgrnorth_paper/planefit_results.txt", ",")
    planes = []
    for i in range(6):
        planes.append(data[i,:])
    #for plane in planes:  print plane
    dots = sc.zeros((6,6))
    for i in range(6):
        for j in range(6):
            dots[i,j] = (sc.arccos(planes[i][0]*planes[j][0] + planes[i][1]*planes[j][1]
                                   + planes[i][2]*planes[j][2]))*deg
    print dots
    print
    errors = sc.zeros((6,6))
    for i in range(6):
        for j in range(6):
            dotty = planes[i][0]*planes[j][0]+planes[i][1]*planes[j][1]+planes[i][2]*planes[j][2]
            factor = -1.0 / sc.sqrt(1 - dotty*dotty)
            print factor
            errors[i,j] = factor*(planes[j][0]*planes[i][4] + planes[i][0]*planes[j][4] +
                                  planes[j][1]*planes[i][5] + planes[i][1]*planes[j][5] +
                                  planes[j][2]*planes[i][6] + planes[i][2]*planes[j][6])*deg
            #if i==3 and j==5:  print dotty
    print errors
    Sgr = [14.655645800014774, 2.2704309216189817, -5.8873772610912614]
    print "# - dist to Sgr Core:"
    for plane in planes:
        print (Sgr[0]*plane[0] + Sgr[1]*plane[1] + Sgr[2]*plane[2] + plane[3]), \
        (Sgr[0]*plane[4] + Sgr[1]*plane[5] + Sgr[2]*plane[6] + plane[7])

def sgr_lam_split(stripes, params):
    # get lambda for stripes
    holder = []
    for i in range(len(stripes)):
        mu, r, theta, phi, wedge = eval(params[i][4]), eval(params[i][5]), \
            eval(params[i][6]), eval(params[i][7]), stripes[i]
        x,y,z = coor.stream2xyz(0.0, 0.0, 0.0, mu, r, theta*deg, phi*deg, wedge)
        Xs, Ys, Zs, lam, beta, r = sl.law_xyz2sgr_sun(x,y,z)
        #if lam > 180.0: lam = lam - 360.0
        #lam = -lam
        print lam, beta, r, stripes[i]
        holder.append([lam, beta, r])
    Xs, Ys, Zs, lam, beta, r = sl.law_xyz2sgr_sun(14.655646, 2.270431, -5.887377)
    print lam, beta, r
    # get lambda range for northern Sgr  15 N stripes, 0-14 -> 9-23
    low_lam, high_lam = 200.0, 314.8  # from Sgr stripes: (200.025901082, 314.707561185)
    size = (high_lam - low_lam) / 15.0
    print "# - LAMBDA BIN SIZE: {0}".format(size)
    for h in range(15):
        array = []
        lam_min, lam_max = low_lam+(size*h), low_lam+(size*(h+1))
        for i in range(15):
            if stripes[i] == 9:  num = "09"
            else:  num = str(stripes[i])
            name = "./sim_streams/sim_stripe_"+num+".txt"
            Sgr = open(name, "r")
            for line in Sgr:
                if line.strip() == "":  continue
                if line.strip()[0] == "#":  continue
                l,b,r,f = line.strip().split()
                x,y,z = coor.lbr2xyz(float(l), float(b), float(r))
                Xs, Ys, Zs, lam, beta, r = sl.law_xyz2sgr_sun(x,y,z)        
                #if lam > high_lam:  high_lam = lam
                #if lam < low_lam:  low_lam = lam
                if lam > lam_min:
                    if lam <= lam_max:
                        array.append([lam, beta, r, eval(f)])
            Sgr.close()
            print "# - Stripe "+num+" successfully read in"
        fi.write_data(sc.array(array), "lamsim_out_"+str(h)+".txt",
                      header="Lambda {0}:{1}, lam, beta, r".format(lam_min, lam_max))
    #print low_lam, high_lam
    # bin stripes in lambda, g
    # correct stripes for detection, FTO efficiencies
    # reflect stars about midpoint to correct for missing data on far end
    print "### --- Done"

def sgr_density(stripes, params):
        # get lambda for stripes
    holder = []
    for i in range(len(stripes)):
        mu, r, theta, phi, wedge = eval(params[i][4]), eval(params[i][5]), \
            eval(params[i][6]), eval(params[i][7]), stripes[i]
        x,y,z = coor.stream2xyz(0.0, 0.0, 0.0, mu, r, theta, phi, wedge)
        Xs, Ys, Zs, lam, beta, r = sl.law_xyz2sgr_sun(x,y,z)
        #if lam > 180.0: lam = lam - 360.0
        #lam = -lam
        print lam, beta, r, stripes[i]
        holder.append([lam, beta, r])
    Xs, Ys, Zs, lam, beta, r = sl.law_xyz2sgr_sun(14.655646, 2.270431, -5.887377)
    print lam, beta, r
    # Get distances as a function of Lamda, using straight pipes to connect points.
    LAB = [] # lambda, A, B such that:  R = A(lam) + B;  C,D st beta = C(lam) + D
    for i in range(15-1):  # all but end
        A = (holder[i+1][2]-holder[i][2])/(holder[i+1][0]-holder[i][0]) # R2-R1/L2-L1
        B = holder[i][2] - A*holder[i][0]  # R1 - A*L1
        C = (holder[i+1][1]-holder[i][1])/(holder[i+1][0]-holder[i][0]) # B2-B1/L2-L1
        D = holder[i][1] - C*holder[i][0]  # B1 - A*L1
        LAB.append([holder[i+1][0], A, B, C, D])
    # get lambda range for northern Sgr  15 N stripes, 0-14 -> 9-23
    low_lam, high_lam = 200.0, 314.8  # from Sgr stripes: (200.025901082, 314.707561185)
    size = (high_lam - low_lam) / 15.0
    print "# - LAMBDA BIN SIZE: {0}".format(size)
    out = []
    for i in range(15):
        killed = 0
        print "# - Lambda wedge {0}:".format(i)
        num = str(i)
        name = "./Lambda_outs/lamsim_out_"+num+".txt"
        data = fi.read_data(name, echo=1)
        count, in_1, in_2, in_0 = 0, 0.0, 0.0, 0.0
        sigma = params[i][8]
        for j in range(len(data[:,0])):
            lam, beta, r = data[j,0], data[j,1], data[j,2]
            count = count + 1
            index = 0
            while data[j,0] > LAB[index][0]:
                index = index + 1
                if (index == len(LAB)-1):  break
            # correct for missing far edge by counting all near edge and doubling
            A, B = LAB[index][1], LAB[index][2]
            if r > (A*lam + B)/sc.cos(beta):  killed=killed+1; continue
            #D = (0.295*lam - 42.601)/sc.cos(beta)
            #if r > D:  killed=killed+1; continue
            # correct for missing side edge by counting inner and doubling
            #C, D = LAB[index][3], LAB[index][4]
            #if beta > (C*lam + D):  continue
            # Correct for completeness
            Deff = 1.0/SDSS_eff(r)
            in_0 = in_0 + 1.0
            in_1 = in_1 + Deff
            in_2 = in_2 + Deff/TO_eff(r)
        out.append([(low_lam+(size*i)+(size*0.5)), count, in_1, in_2, in_2/size])
        print "###   Number of stars in stream:  {0}".format(count)
        print "###   KILLED: {0}".format(killed)
        print "###   Remaining Stars:    {0}".format(in_0*1.0)
        print "###   Corrected for SDSS Eff.:    {0}".format(in_1*2.0)
        print "###   Corrected for SDSS/FTO Eff: {0}".format(in_2*2.0)
    #print out
    return sc.array(out)

def sgr_density_south(stripes, params):
    out = []
    for i in range(15, 18):
        mu, r, theta, phi, wedge = eval(params[i][4]), eval(params[i][5]), \
            eval(params[i][6]), eval(params[i][7]), stripes[i]
        print mu, r, theta, phi
        x,y,z = coor.stream2xyz(0.0, 0.0, 0.0, mu, r, theta, phi, wedge)
        Xs, Ys, Zs, lam, beta, r = sl.law_xyz2sgr_sun(x,y,z)
        print lam, beta, r, stripes[i]
        name = "./sep_lbr/s1-"+str(stripes[i])+".txt"
        data = fi.read_data(name, echo=1)
        wedge_norm = coor.stripe_normal(stripes[i])
        plane_norm = [0.03176573688051381, 0.98687684085235916, 0.15831942063343135]
        stream_norm = [sc.cos(phi)*sc.sin(theta), sc.sin(phi)*sc.sin(theta), sc.cos(theta)]
        print wedge_norm, stream_norm
        dotty = wedge_norm[0]*stream_norm[0] + wedge_norm[1]*stream_norm[1] + wedge_norm[2]*stream_norm[2]
        #dotty2 = plane_norm[0]*wedge_norm[0] + plane_norm[1]*stream_norm[1] + plane_norm[2]*stream_norm[2]
        if dotty > 1.0:  print "!!!  Not a valid dot-product: {0}".format(dotty)
        H = abs(2.5 / dotty)
        #H2 = 2.5 / sc.sqrt(1 - (dotty2*dotty2))
        print "H = ", H, dotty, sc.arccos(dotty)*180.0/ma.pi #, H2
        count, in_1, in_2 = 0, 0.0, 0.0
        lam_low, lam_high = 360.0, 0.0
        for j in range(len(data[:,0])):
            l, b, r = data[j,0], data[j,1], data[j,2]
            count = count + 1
            # Correct for completeness
            Deff = 1.0/SDSS_eff(r)
            in_1 = in_1 + Deff
            in_2 = in_2 + Deff + 1.0/TO_eff(r)
            x,y,z = coor.lbr2xyz(l,b,r)
            Xs, Ys, Zs, lam2, beta2, r2 = sl.law_xyz2sgr_sun(x,y,z)
            if lam2 < lam_low:  lam_low = lam2
            if lam2 > lam_high:  lam_high = lam2
        out.append([lam, count, in_1, in_2, in_2/H])
        print "--- Lambda Range: {0}:{1}".format(lam_low, lam_high)
        print "###   Number of stars in stream:  {0}".format(count)
        print "###   Corrected for SDSS Eff.:    {0}".format(in_1)
        print "###   Corrected for SDSS/FTO Eff: {0}".format(in_2)
    for o in out:  print o
    return sc.array(out)

def sgr_sim_density(stripes, params):
    SGA = fi.read_data("../sgrnorth_paper/stream_gen_analysis.txt", ",")
    low_lam, high_lam = 200.0, 314.8  # from Sgr stripes: (200.025901082, 314.707561185)
    size = (high_lam - low_lam) / 15.0
    N_stars = [19386, 18734, 11096, 17044, 18736, 15409, 12519, 12248, 8853, 7328,
               5479, 4450, 3486, 2425, 971, 9511, 16119, 16603]
    south = [79, 82, 86]
    out = []
    for i in range(15):
        print "# - Lambda wedge {0}:".format(i)
        num = str(i)
        name = "./sim_streams/lamsim_out_"+num+".txt"
        data = fi.read_data(name, echo=1)
        total = len(data[:,0])
        other = sc.sum(data[:,3])
        length = 7.653333333  # bin size
        out.append([(low_lam+(size*i)+(size*0.5)), total-other, total, 0.0, total/length])
    for i in range(15,18):
        print "# - Southern wedge {0}:".format(south[i-15])
        data = fi.read_data("./sim_streams/sim_stripe_"+str(south[i-15])+".txt")
        total = len(data[:,0])
        other = sc.sum(data[:,3])
        length = SGA[i,3]
        mu, r, theta, phi, wedge = eval(params[i][4]), eval(params[i][5]), \
            eval(params[i][6]), eval(params[i][7]), stripes[i]
        x,y,z = coor.stream2xyz(0.0, 0.0, 0.0, mu, r, theta, phi, wedge)
        Xs, Ys, Zs, lam, beta, r = sl.law_xyz2sgr_sun(x,y,z)
        out.append([lam, total-other, total, 0.0, total/length])
    fi.write_data(sc.array(out), "sgrNorth_density_sim.txt", ",", " lam(mid), count, SDSS Deff, SDSS/FTO Deff, per degree")


def SDSS_eff(r):
    g0 = coor.getg(r, 4.2)
    bottom = sc.exp(1.6171*(g0-23.5877)) + 1.0
    return 0.9402/bottom

def do_SDSS_eff():
    x = sc.arange(0.0, 50.0, 1.0)
    y = TO_eff(x)
    #y = SDSS_eff(x)
    plt.figure()
    plt.plot(x,y)
    plt.show()

# Initialize params for turnoff efficiency
ay = [8.55878159e+00, -1.04891551e+01, 3.51630757e+00, -2.29741062e-01, 6.72278105e-03, -1.01910181e-04, 7.82787167e-07, -2.41452056e-09]
ar = [5.61945007e+02, -1.67343282e+01, 1.09325822e-01, 1.34993610e-03, -1.42044161e-05]
af = ay
for i in range(len(ar)):
    af[i] = ay[i] + ar[i]
def Nth_order_polynomial(x, params):
    #A polynomial of order=number of parameters - 1
    y = params[0]
    for N in range(1,len(params)):
        holder = 1.0
        for j in range(N):
            holder = holder*x
        y = y + params[N]*holder
    return y/532.0
def TO_eff(r):
    return Nth_order_polynomial(r, af)

def show_lambda():
    data = fi.read_data("./sep_lbr/s1-79.txt")
    x,y,z = coor.lbr2xyz(data[:,0],data[:,1],data[:,2])
    Xs,Ys,Zs,lam1,beta,r = sl.law_xyz2sgr_sun(x,y,z)
    data = fi.read_data("./sep_lbr/s1-82.txt")
    x,y,z = coor.lbr2xyz(data[:,0],data[:,1],data[:,2])
    Xs,Ys,Zs,lam2,beta,r = sl.law_xyz2sgr_sun(x,y,z)
    data = fi.read_data("./sep_lbr/s1-86.txt")
    x,y,z = coor.lbr2xyz(data[:,0],data[:,1],data[:,2])
    Xs,Ys,Zs,lam3,beta,r = sl.law_xyz2sgr_sun(x,y,z)
    plt.figure()
    plt.hist(lam3, 20)
    plt.hist(lam2, 20)
    plt.hist(lam1, 20)
    plt.show()
    plt.close('all')

def width_plot(stripes, params):
    data = fi.read_data("../sgrnorth_paper/sgrNorth_density_sim.txt", ",")
    #Ben = fi.read_data("../sgrnorth_paper/bensgr_lam_bins.dat")
    #Law = fi.read_data("../sgrnorth_paper/lawsgr_lam_bins.dat")
    Law = fi.read_data("../sgrnorth_paper/lawsgrTRX_lam_binsC.dat")
    sph = fi.read_data("../sgrnorth_paper/lawsgrSPH_lam_binsC.dat")
    obl = fi.read_data("../sgrnorth_paper/lawsgrOBL_lam_binsC.dat")
    pro = fi.read_data("../sgrnorth_paper/lawsgrPRO_lam_binsC.dat")
    LM10 = fi.read_data("../sgrnorth_paper/lawsgrTRX_widthsC.dat")  #2
    wp = fi.read_data("../sgrnorth_paper/lawsgrPRO_widthsC.dat")
    ws = fi.read_data("../sgrnorth_paper/lawsgrSPH_widthsC.dat")
    wo = fi.read_data("../sgrnorth_paper/lawsgrOBL_widthsC.dat")
    errors = sc.array([0.775105, 0.47734, 0.836774, 1.882841, 0.932078, 1.531246,
                       1.653404, 0.593513, 0.853309, 0.469178, 3.053077, 0.314356,
                       0.388961, 1.21, 2.330186, 0.142532, 0.989839, 0.541965])
    sigma, slam = [], []
    for i in range(len(stripes)):
        mu, r, theta, phi, wedge = eval(params[i][4]), eval(params[i][5]), \
            eval(params[i][6]), eval(params[i][7]), stripes[i]
        #print mu, r, theta, phi
        x,y,z = coor.stream2xyz(0.0, 0.0, 0.0, mu, r, theta, phi, wedge)
        Xs, Ys, Zs, lam, beta, r = sl.law_xyz2sgr_sun(x,y,z)
        slam.append(lam)
    for param in params:
        sigma.append(2.35482*eval(param[8]))
    fig = plt.figure(1)
    plt.subplots_adjust(hspace=0.001, wspace=0.001)
    # density along stream plot
    sp1 = plt.subplot(211)
    color = ['w', 'w', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'w', 'w', 'w', 'k', 'k', 'k']
    #plt.scatter(data[:,0], data[:,4]*1.0, facecolor=color, marker="s",
    #            label="Expected Counts")
    plt.scatter(data[:,0], data[:,4], facecolor=color, marker="^", label="__nolabel__")
    plt.scatter(data[:,0], data[:,1]*(data[:,4]/data[:,2]), facecolor=color,
                marker="o", label="__nolabel__")
    #plt.errorbar([315.0, 315.0, 315.0, 312.0, 312.0, 312.0], [1000, 2000, 3000, 1000, 2000, 3000],
    #            yerr=[39.2, 51.5, 80.0, 33.2, 43.6, 53.8], fmt="o", markersize=1,
    #            mfc="red", ecolor="red", label="Typical Errors") #71.8
    #plt.plot(Ben[:,0], Ben[:,1]*0.2, "k--") #0.1
    a, b = 21, 23
    plt.plot(obl[:a,0], obl[:a,1]*7.0, "r:", label="Oblate")
    plt.plot(obl[b:,0], obl[b:,1]*7.0, "r:", label=None)
    plt.plot(sph[:a,0], sph[:a,1]*7.0, "k-.", label="Spherical", linewidth=1.5)
    plt.plot(sph[b:,0], sph[b:,1]*7.0, "k-.", label=None, linewidth=1.5)
    plt.plot(pro[:a,0], pro[:a,1]*7.0, "b--", label="Prolate")
    plt.plot(pro[b:,0], pro[b:,1]*7.0, "b--", label=None)
    plt.plot(Law[:a,0], Law[:a,1]*6.0, "k-", label="L&M10 Nbody") # 1.5
    plt.plot(Law[b:,0], Law[b:,1]*6.0, "k-", label=None) # 1.5
    #sp1.annotate('Triaxial', xy=(178, 3200), xytext=(210, 3650), fontsize=8,
    #             arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=3)  )
    #sp1.annotate('Prolate', xy=(192, 3000), xytext=(215, 3350), fontsize=8,
    #             arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=3)  )
    #sp1.annotate('Spherical', xy=(199, 2750), xytext=(220, 3050), fontsize=8,
    #             arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=3)  )
    #sp1.annotate('Oblate', xy=(205, 2500), xytext=(225, 2800), fontsize=8,
    #             arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=3)  )
    #plt.scatter(Ben[:,0], Ben[:,1]*0.5, facecolor="blue", marker="x")
    #plt.scatter(Law[:,0], Law[:,1], facecolor="red", marker="+")
    spoof1 = plt.scatter(180.0, -10.0, facecolor='k', marker="^", label="Expected Counts") #spoof
    spoof2 = plt.scatter(180.0, -10.0, facecolor='k', marker="o", label="Raw Counts") #spoof
    plt.setp(sp1.get_xticklabels(), visible=False)
    plt.ylim(0.0, 6000.0)
    plt.ylabel("Stars per degree", fontsize=12)
    # legend
    leg1 = sp1.legend(loc='upper center', numpoints=2)
    for t in leg1.get_texts():
        t.set_fontsize(10)
    # ----------------  FWHM v lambda plot --------------------
    sp2 = plt.subplot(212, sharex=sp1)
    g_mean = np.sqrt(LM10[:,3]*LM10[:,4])  #geometric mean
    plt.errorbar(slam, sigma, yerr=2.35482*errors, fmt="d", mfc="black", ecolor="black", label="MLE Fits")
    #plt.scatter(slam, sigma, facecolor="black", marker="o", label=None)
    plt.plot(LM10[:21,0], 2.35482*g_mean[:21], "k-", label=None)
    plt.plot(LM10[21:,0], 2.35482*g_mean[21:], "k-", label=None)
    # NEW STUFF
    aa, bb = 21, 21
    plt.plot(wp[:aa,0], 2.35482*np.sqrt(wp[:aa,3]*wp[:aa,4]), "b--", label=None)
    plt.plot(wp[bb:,0], 2.35482*np.sqrt(wp[bb:,3]*wp[bb:,4]), "b--", label=None)
    plt.plot(ws[:aa,0], 2.35482*np.sqrt(ws[:aa,3]*ws[:aa,4]), "k-.", label=None)
    plt.plot(ws[bb:,0], 2.35482*np.sqrt(ws[bb:,3]*ws[bb:,4]), "k-.", label=None)
    plt.plot(wo[:aa,0], 2.35482*np.sqrt(wo[:aa,3]*wo[:aa,4]), "r:", label=None)
    plt.plot(wo[bb+3:,0], 2.35482*np.sqrt(wo[bb+3:,3]*wo[bb+3:,4]), "r:", label=None)
    #plt.scatter(LM10[:,0], 2.35482*g_mean, facecolor='black', marker='+', label="L&M10 Nbody")
    # legend
    leg2 = sp2.legend(loc='lower right', numpoints=1)
    for t in leg2.get_texts():
        t.set_fontsize(8)
    plt.xlabel(r"$\Lambda_{\odot}$")
    plt.ylabel(r"FWHM (kpc)")
    plt.ylim(0.0, 20.0)
    plt.xlim(320.0, 70.0) # 70, 320
    #plt.xticks(sc.arange(80.0, 320.0, 20.0))
    plt.show()

def compare_MW():
    stripes1, params1 = load_params("../sgrnorth_paper/results_BG_easyread.txt")
    stripes2, params2 = load_params("../sgrnorth_paper/results_MW_easyread.txt")
    errors=fi.read_data("../sgrnorth_paper/results_errors_final.txt")
    table = []
    for i in range(15):
        if float(params2[i][0]) > float(params1[i][0]):  line = "YES "
        else:  line = "NO  "
        table_l = str(stripes1[i])+" & "    
        for j in range(1,9):
            #print params1[i][j], params2[i][j], errors[i,(j-1)]
            diff = abs(float(params1[i][j])-float(params2[i][j])) / errors[i,(j-1)]
            #diff = 100.0*abs(float(params1[i][j])-float(params2[i][j])) / float(params1[i][j])
            line = line + str(round(diff,2)) + "\t"  #FOR JUST Dsigmas
            param = float(params2[i][j])
            table_l = table_l + str(round(param,2)) + "/" + str(round(diff,2)) + " & "
        table.append(table_l[:-2])
        print stripes1[i], stripes2[i], line
    for item in table:
        print item

def make_lbr_MRT():
    import glob
    out = []
    stripes = ["09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
               "20", "21", "22", "23", "79", "82", "86"]
    backs = glob.glob("./sep_lbr/Hold/*.txt")
    #for stripe in stripes:
        #infile = open("./sep_lbr/s1-"+stripe+".txt", "r")
    for back in backs:
        infile = open(back, "r")
        for line in infile:
            l = "  "
            holder = line.strip().split()
            l = l + holder[0]
            l = l + "\t" + holder[1]
            g = coor.getg(float(holder[2]))
            g = str(round(g,2))
            if len(g) < 5:
                g = g + "0"
            l = l + "\t" + g 
            out.append(l)
        infile.close()
    #outfile = open("sgr_separated_stars_MRT.txt", "w")
    outfile = open("back_separated_stars_MRT.txt", "w")
    for o in out:
        outfile.write(o+"\n")
    outfile.close()

def compare_halos():
    A, B = 1.0, 1.0
    x,y = 8.5, 0.0
    z = np.arange(0.0, 30.1, 1.0)
    # power law (Juric+08)
    n_h, q1 = 2.8, 0.64
    power_law = A*(8.5 / (np.sqrt(x*x + y*y) + ( (z*z)/(q1*q1) ) ) )**n_h
    # Hernquist
    r0, q2 = 6.73, 0.53
    rr = np.sqrt(x*x + y*y + ((z*z)/(q2*q2)) )
    denom = rr*(rr + r0)*(rr + r0)*(rr + r0)
    B = np.ma.min(denom)
    hernquist = B / denom
    plt.figure(1)
    plt.plot(z, power_law, "k-")
    plt.plot(z, hernquist, "k--")
    plt.show()

def koposov_geti():
    """ get i values from Koposov 2012, updated with Table 2 errata """
    lams = [92.5, 97.5, 102.5, 107.5, 112.5, 117.5, 122.5, 127.5]
    mags_old = [16.69, 16.72, 16.86, 17.01, 17.17, 17.31, 17.27, 17.28]
    mags_new = [17.04, 17.07, 17.21, 17.36, 17.52, 17.66, 17.62, 17.63]
    slope = 0.023
    print "Original Table:"
    for mags in (mags_old, mags_new):
        # For lambda = 75.0, 85.0, get line from gradient
        m = slope
        b = mags[0] + 0.6 - (lams[0]*m)
        #m = (16.72-16.69)/(97.5-92.5)
        #b = 16.69+0.6 - (92.5*m)
        i = 75.0*m + b
        print "lam=75.0, i={0}, m={1}, b={2}".format(i,m,b)
        i = 85.0*m + b
        print "lam=85.0, i={0}, m={1}, b={2}".format(i,m,b)
        # For lambda = 95.0, get line from adjacent points
        m = (mags[1]-mags[0])/(lams[1]-lams[0])
        b = mags[1]+0.6 - (lams[1]*m)
        i = 95.0*m + b
        print "lam=95.0, i={0}, m={1}, b={2}".format(i,m,b)
        # For lambda = 105.0, get line from adjacent points
        m = (mags[3]-mags[2])/(lams[3]-lams[2])
        b = mags[3]+0.6 - (lams[3]*m)
        i = 105.0*m + b
        print "lam=115.0, i={0}, m={1}, b={2}".format(i,m,b)
        # For lambda = 115.0, get line from adjacent points
        m = (mags[5]-mags[4])/(lams[5]-lams[4])
        b = mags[5]+0.6 - (lams[5]*m)
        i = 115.0*m + b
        print "lam=115.0, i={0}, m={1}, b={2}".format(i,m,b)
        print "\n Errata Table:"
    

if __name__ == "__main__":
    stripes, params = load_params("../sgrnorth_paper/results_BG_easyread.txt")
    #results_to_table(stripes, params)
    #sgr_xyz_plots(stripes, params)
    width_plot(stripes, params)
    #get_errors()
    #get_distances()
    #make_coord_table(stripes, params)
    #astro_test()
    #dotties()
    #plane_fit_table()
    #sgr_lam_split(stripes, params)
    #sgr_sim_density(stripes, params)
    #out = sgr_density(stripes, params)
    #print out
    #fi.write_data(out, "sgrNorthsim_density.txt", ",", " lam(mid), count, SDSS Deff, SDSS/FTO Deff, per degree")
    #sgr_density_south(stripes, params)
    #show_lambda()
    #show_vecs()
    #do_SDSS_eff()
    #sda.integrate_stripes(stripes, params)
    #make_table_werrors(stripes, params)
    #compare_MW()
    #make_lbr_MRT()
    #compare_halos()
    #koposov_geti()
