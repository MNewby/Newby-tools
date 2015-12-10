import math as ma
import numpy as np
import scipy as sc

""" Law Sgr code converted into Python
    Matthew Newby, RPI, July 5, 2012 """


# intialize constants and matrices
dsun = 8.0  #Formerly 8.5 - Oct 1, 2012
radpdeg=3.141592653589793/180.
#// Define the Euler angles
phi=(180.0+3.75)*radpdeg
theta=(90.0-13.46)*radpdeg
psi=(180.0+14.111534)*radpdeg
#// Define the rotation matrix from the Euler angles
rot11=sc.cos(psi)*sc.cos(phi) - sc.cos(theta)*sc.sin(phi)*sc.sin(psi)
rot12=sc.cos(psi)*sc.sin(phi) + sc.cos(theta)*sc.cos(phi)*sc.sin(psi)
rot13=sc.sin(psi)*sc.sin(theta)
rot21=-sc.sin(psi)*sc.cos(phi) - sc.cos(theta)*sc.sin(phi)*sc.cos(psi)
rot22=-sc.sin(psi)*sc.sin(phi) + sc.cos(theta)*sc.cos(phi)*sc.cos(psi)
rot23=sc.cos(psi)*sc.sin(theta)
rot31=sc.sin(theta)*sc.sin(phi)
rot32=-sc.sin(theta)*sc.cos(phi)
rot33=sc.cos(theta)
#Define rotation matrices
rotM = np.matrix([[rot11, rot12, rot13],[rot21, rot22, rot23],[rot31,rot32,rot33]])
rotI = rotM.I  #inverse of M

def law_xyz2sgr_sun(X,Y,Z, Xsun=dsun):
    """ From the Law sgr transformation code:  http://www.astro.virginia.edu/~srm4n/Sgr/code.html
    // Modified to move contants outside of function
    // Transform positions from standard right handed Galactocentric XYZ to
    // the heliocentric Sgr system (lambda=0 at Sgr)
    // Input must be in kpc of the form X Y Z
    // Output is in kpc and degrees, of the form X_Sgr Y_Sgr Z_Sgr r lambda beta"""
    X = X + Xsun #// Transform the input system to heliocentic right handed coordinates
    #// Calculate X,Y,Z,distance in the Sgr system
    galXYZ = np.matrix([X,Y,Z])
    sgrXYZ = rotM*galXYZ.T
    r = np.linalg.norm(sgrXYZ)
    Xs, Ys, Zs = sgrXYZ[0,0], sgrXYZ[1,0], -sgrXYZ[2,0]  # Z is flipped due to beta pole orientation
    #// Calculate the angular coordinates lambda,beta
    lam = sc.arctan2(Ys,Xs)/radpdeg
    if type(lam) != type(sc.zeros(0)):
        if (lam < 0):  lam = lam + 360.0
    else:
        for i in range(len(lam)):
            if (lam[i] < 0):  lam[i] = lam[i] + 360.0
    beta = sc.arcsin(Zs/sc.sqrt(Xs*Xs + Ys*Ys + Zs*Zs))/radpdeg
    return Xs, Ys, Zs, lam, beta, r

def law_sgr2xyz_sun(lam, beta, rr, Xsun=dsun):
    """ inverse of previous function """
    Xs = rr*sc.cos(lam*radpdeg)*sc.cos(beta*radpdeg)
    Ys = rr*sc.sin(lam*radpdeg)*sc.cos(beta*radpdeg)
    Zs = -1.0*rr*sc.sin(beta*radpdeg)  #beta pole is flipped
    sgrXYZ = np.matrix([Xs, Ys, Zs])
    galXYZ = rotI*sgrXYZ.T
    return galXYZ[0,0]-Xsun, galXYZ[1,0], galXYZ[2,0]

def law_xyz2sgr_sun_old(X,Y,Z, Xsun=dsun):
    """ From the Law sgr transformation code:  http://www.astro.virginia.edu/~srm4n/Sgr/code.html
    // Transform positions from standard right handed Galactocentric XYZ to
    // the heliocentric Sgr system (lambda=0 at Sgr)
    // Input must be in kpc of the form X Y Z
    // Output is in kpc and degrees, of the form X_Sgr Y_Sgr Z_Sgr r lambda beta"""
    radpdeg=3.141592653589793/180.
    #// Define the Euler angles
    phi=(180.0+3.75)*radpdeg
    theta=(90.0-13.46)*radpdeg
    psi=(180.0+14.111534)*radpdeg
    #// Define the rotation matrix from the Euler angles
    rot11=sc.cos(psi)*sc.cos(phi) - sc.cos(theta)*sc.sin(phi)*sc.sin(psi)
    rot12=sc.cos(psi)*sc.sin(phi) + sc.cos(theta)*sc.cos(phi)*sc.sin(psi)
    rot13=sc.sin(psi)*sc.sin(theta)
    rot21=-sc.sin(psi)*sc.cos(phi) - sc.cos(theta)*sc.sin(phi)*sc.cos(psi)
    rot22=-sc.sin(psi)*sc.sin(phi) + sc.cos(theta)*sc.cos(phi)*sc.cos(psi)
    rot23=sc.cos(psi)*sc.sin(theta)
    rot31=sc.sin(theta)*sc.sin(phi)
    rot32=-sc.sin(theta)*sc.cos(phi)
    rot33=sc.cos(theta)
    # X=-X; // Make the input system right-handed
    X = X + Xsun #// Transform the input system to heliocentic right handed coordinates
    #// Calculate X,Y,Z,distance in the Sgr system
    Xs=rot11*X + rot12*Y + rot13*Z
    Ys=rot21*X + rot22*Y + rot23*Z
    Zs=rot31*X + rot32*Y + rot33*Z
    r=sc.sqrt(Xs*Xs + Ys*Ys + Zs*Zs)
    Zs=-Zs
    #// Calculate the angular coordinates lambda,beta
    lam = sc.arctan2(Ys,Xs)/radpdeg
    if type(lam) != type(sc.zeros(0)):
        if (lam < 0):  lam = lam + 360.0
    else:
        for i in range(len(lam)):
            if (lam[i] < 0):  lam[i] = lam[i] + 360.0
    beta = sc.arcsin(Zs/sc.sqrt(Xs*Xs + Ys*Ys + Zs*Zs))/radpdeg
    return Xs, Ys, Zs, lam, beta, r


def law_xyz2sgr_gal(X,Y,Z, Xsun=dsun):
    """ From the Law sgr transformation code:  http://www.astro.virginia.edu/~srm4n/Sgr/code.html
    // Transform positions from standard left handed Galactocentric XYZ to  [RIGHT-HANDED]
    // the Galactocentric Sgr system (lambda=0 at the Galactic plane)
    // Input must be in kpc of the form X Y Z
    // Output is in kpc and degrees, of the form X_Sgr,GC Y_Sgr,GC Z_Sgr,GC d_GC lambda_GC beta_GC
    // Note that d is distance from Galactic Center"""
    radpdeg = 3.141592653589793/180.0
    #// Define the Euler angles
    phi = (180+3.75)*radpdeg
    theta = (90-13.46)*radpdeg
    psiGC = (180+21.604399)*radpdeg
    #// Rotation angle of phiGC past 180degrees is a useful number
    ang = 21.604399*radpdeg
    #// Note that the plane does not actually include the G.C., although it is close
    xcenter = -8.5227
    ycenter = -0.3460
    zcenter = -0.0828
    #// Define the rotation matrix from the Euler angles
    GCrot11 = sc.cos(psiGC)*sc.cos(phi)-sc.cos(theta)*sc.sin(phi)*sc.sin(psiGC)
    GCrot12 = sc.cos(psiGC)*sc.sin(phi)+sc.cos(theta)*sc.cos(phi)*sc.sin(psiGC)
    GCrot13 = sc.sin(psiGC)*sc.sin(theta)
    GCrot21 = -sc.sin(psiGC)*sc.cos(phi)-sc.cos(theta)*sc.sin(phi)*sc.cos(psiGC)
    GCrot22 = -sc.sin(psiGC)*sc.sin(phi)+sc.cos(theta)*sc.cos(phi)*sc.cos(psiGC)
    GCrot23 = sc.cos(psiGC)*sc.sin(theta)
    GCrot31 = sc.sin(theta)*sc.sin(phi)
    GCrot32 = -sc.sin(theta)*sc.cos(phi)
    GCrot33 = sc.cos(theta)
    #X = -X  #// Make the input system right-handed  ALREADY IS
    X = X + Xsun  #// Transform the input system to heliocentric right handed coordinates
    #// Calculate Z,distance in the SgrGC system
    Temp=GCrot11*(X+xcenter)+GCrot12*(Y-ycenter)+GCrot13*(Z-zcenter)
    Temp2=GCrot21*(X+xcenter)+GCrot22*(Y-ycenter)+GCrot23*(Z-zcenter)
    Zs=GCrot31*(X+xcenter)+GCrot32*(Y-ycenter)+GCrot33*(Z-zcenter)
    d=sc.sqrt(Temp*Temp+Temp2*Temp2+Zs*Zs)
    Zs=-Zs;
    #// Calculate the angular coordinates lambdaGC,betaGC
    Temp3 = sc.arctan2(Temp2,Temp)/radpdeg;
    if (Temp3 < 0.0): Temp3 = Temp3 + 360.0
    Temp3 = Temp3 + ang/radpdeg
    if (Temp3 > 360.0): Temp3 = Temp3 - 360.0
    lam = Temp3
    beta = sc.arcsin(Zs/sc.sqrt(Temp*Temp+Temp2*Temp2+Zs*Zs))/radpdeg
    #// Calculate X,Y in the SgrGC system
    Xs=Temp*sc.cos(ang) - Temp2*sc.sin(ang)
    Ys=Temp*sc.sin(ang) + Temp2*sc.cos(ang)
    r = sc.sqrt(Xs*Xs + Ys*Ys + Zs+Zs)
    return Xs, Ys, Zs, lam, beta, r

"""
if __name__ == "__main__":
    X, Y, Z = -8.0, 0.0, 20.0
    aa = law_xyz2sgr_sun_old(X,Y,Z)
    bb = law_xyz2sgr_sun(X,Y,Z)
    print aa
    print bb
    if aa == bb:  print "YES!"
    print law_sgr2xyz_sun(aa[3], aa[4], aa[5])
"""
