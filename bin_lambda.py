import numpy as np

size = 7.653333333
reference = 200.0
num_bins = 1 + (360.0 / size)
rad = np.pi/180.0

centers = np.arange((size/2.0), 367.5, size)
counts = np.zeros(num_bins, int)
degrees = np.zeros(361, int)
slices = []
for i in range(len(centers)):
    slices.append([])

DO = "OBL"  # TRX, PRO, SPH, or OBL

# initialize lambda, W cuts, filenames, Pcol limit, Pcol and R indices.
# Files from Law's Sgr simulations (Law+05, Law+10)
if DO == "PRO":
    p1 = [0.0, 128.0]
    p2 = [185.0, -200.0]
    p3 = [316.0, 220.0]
    filename = "prolate.dat"
    Pcol, Pi, ri = 2.0, 15, 13
if DO == "SPH":
    p1 = [0.0, 160.0]
    p2 = [170.0, -220.0]
    p3 = [320.0, 280.0]
    filename = "spherical.dat"
    Pcol, Pi, ri = 2.0, 15, 13
if DO == "OBL":
    p1 = [0.0, 160.0]
    p2 = [160.0, -220.0]
    p3 = [260.0, 60.0]
    filename = "oblate.dat"
    Pcol, Pi, ri = 2.0, 15, 13
if DO == "TRX":
    filename = "SgrTriax_DYN.dat"
    Pcol, Pi, ri = 3.0, 24, 18

# Initialize axisyemmtric cuts;  W = a*lam + b
if DO != "TRX":
    a1 = (p1[1]-p2[1])/(p1[0]-p2[0])
    b1 = p1[1] - a1*p1[0]
    a2 = (p2[1]-p3[1])/(p2[0]-p3[0])
    b2 = p2[1] - a2*p2[0]
    print a1, b1, a2, b2


# Actually do the binning
infile = open(filename, "r")
print "# - File Opened"
for line in infile:
    if line.strip()[0]=="#":  continue
    if line.strip()=="":  continue
    # "3.0" Chosen to include most contiguous particles - lower removes
    #     particles from main stream area, higher adds particles off of the main
    #     stream.  4 seemed optimal, but the additional material was much wider.
    if float(line.split()[Pi]) > Pcol:  continue
    lam, beta, r = float(line.split()[0]), float(line.split()[1]), float(line.split()[ri]) 
    if lam > 160.0 and lam < 185.0:  continue  # Anti-center cut
    """ if Triaxial - uses lmflag and position in sky to separate major streams """
    """ else axisymmetric - uses linear cuts in lambda, W to distiguish streams """
    if DO == "TRX":
        flag = int(line.split()[25])  #lmFlag
        if lam > 180.0:  #North
            if flag < 0:  continue  #If trailing and North, skip
        else:  #South
            if flag > 0:  continue  #If leading and South, skip
    else:
        # Kill below first line
        if lam > p1[0] and lam < p2[0]:  # In Range?
            if float(line.split()[12]) < (a1*lam + b1):  continue
        # kill above second line
        if lam > p2[0] and lam < p3[0]:  # In Range?  
            if float(line.split()[12]) > (a2*lam + b2):  continue
    # If it survived, add star to bins
    index = int(lam/size)
    counts[index] += 1
    index = int(lam)
    degrees[index] += 1
    index = int(lam/size)
    slices[index].append([lam, beta, r])
infile.close()
print "# - File Closed"

writefile = open("lawsgr"+DO+"_lam_binsC.dat", "a")
writefile.write("# bin, centers (lambda), counts from Law/Maj Sgr model, bin size = "+ str(size) + "\n" )
writefile.write("\n")
for i in range(len(centers)):
    writefile.write(str(centers[i])+"\t"+str(counts[i])+"\n")
writefile.close()

writefile = open("lawsgr"+DO+"_deg_binsC.dat", "a")
writefile.write("# bin, centers (lambda), counts from Law/Maj Sgr model, bins size=1.0\n")
writefile.write("\n")
for i in range(len(degrees)):
        writefile.write(str(0.5+1.0*i)+"\t"+str(degrees[i])+"\n")
writefile.close()

writefile = open("lawsgr"+DO+"_widthsC.dat", "a")
writefile.write("# - lambda, x-mean, y-mean, x-stdev, y-stdev, geometric_mean; slice number\n")
writefile.write("\n")
for i in range(len(slices)):
    if slices[i] == []:  writefile.write("# Skipping "+str(i)+"\n"); continue
    stuff = np.array(slices[i])
    X, Y = stuff[:,2]*np.cos(stuff[:,1]*rad), stuff[:,2]*np.sin(stuff[:,1]*rad)
    xc, yc, xs, ys = np.mean(X), np.mean(Y), np.std(X), np.std(Y)
    g_mean = np.sqrt(xs*ys)
    newline = str(centers[i])+"\t"+str(xc)+"\t"+str(yc)+"\t"+str(xs)+"\t"+str(ys)+"\t"+str(g_mean)+"\t"+str(i)+"\n"
    writefile.write(newline)
writefile.close()
