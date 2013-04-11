import math as m

"""This program loads coordinate conversion tools.

Created By Matthew Newby (RPI), 2008
"""

""" !!! DEPRECATED, THIS CODE HAS BEEN MOVED TO 'astro_coordinates.py' and improved"""

print 'Coordinate Tools Loaded:'
print '"getr (g, M=4.2)" converts g (mag) into r (kpc)'
print '"getg (r, M=4.2)" converts r (kpc) into g (mag)'
print '"getM (g,r)" converts g into M, given r'
print '"raDeg (raHr, raMin, raSec, dec=0.0)" converts time RA into degrees'

#converts a magnitude into a distance (kpc)
def getr (g, M=4.2):
    r = ( ( 10.**( (g-M)/5. ) )/ 100. )
    #print 'coverting (magnitudes) g=', g, 'into (kpc) r, using M=', M
    return r 

#converts a distance (kpc) into a magnitude
def getg (r, M=4.2):
    g = M + 5.*(m.log10(r*1000) - 1.)
    #print 'converting (kpc) r=', r, 'into (magnitudes) g, using M=', M
    return g

#converts an apparent magnitude to absolute magnitude, given the distance
def getM (g, d):
    M = g - 5.*(m.log10(d*1000) - 1.)
    #
    return M

#changes ra from time units to degrees, dec in degrees
def raDeg (raHr, raMin, raSec, dec=0.0):
    raTime = raHr + (raMin/60.0) + (raSec/3600.0)
    newra = raTime*(15.0)*(m.cos(m.radians(dec)) )
    #should be dec != 0.0?  Works if dec =0.0...
    #15.0 = (360 degrees)/(24 hours)
    print 'converted ra, in degrees:'
    return newra

#difference in ra (time) in ra (deg)
#(ra1-ra2)*15.0*m.cos(dec)
