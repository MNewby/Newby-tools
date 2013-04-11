import numpy as np
import scipy as sc

rad = np.pi/180.0 

size = 7.653333333
reference = 200.0

centers = np.arange((size/2.0), 367.5, size)
slices = []
for i in range(len(centers)):
	slices.append([])

infile = open("SgrTriax_DYN.dat", "r")
#infile = open("spherical.dat", "r")
for line in infile:
	if line.strip()[0]=="#":  continue
	if line.strip()=="":  continue
	# "3.0" Chosen to include most contiguous particles - lower removes particles from main stream area, higher adds particles off of the main stream.  4 seamed optimal, but the additional material was much wider.
	if float(line.split()[24]) > 4.0:  continue  
	#if float(line.split()[15]) > 1.0:  continue  
	lam, beta, r = float(line.split()[0]), float(line.split()[1]), float(line.split()[18])
	#lam, beta, r = float(line.split()[0]), float(line.split()[1]), float(line.split()[13])
	""" """
	# This part checks to see if the stream is leading or trailing
	flag = int(line.split()[25])
	if lam > 180.0 and lam < 360.0:
		if flag < 0:  continue
	else:
		if flag > 0:  continue
	index = int(lam/size)
	slices[index].append([lam, beta, r])
infile.close()
print "# - Done reading in data"

print "# - lambda, x-mean, y-mean, x-stdev, y-stdev, distance(width); slice number"
for i in range(len(slices)):
	if slices[i] == []:  print "Skipping "+str(i); continue
	stuff = sc.array(slices[i])
	X, Y = stuff[:,2]*sc.cos(stuff[:,1]*rad), stuff[:,2]*sc.sin(stuff[:,1]*rad)
	xc, yc, xs, ys = sc.mean(X), sc.mean(Y), sc.std(X), sc.std(Y)
	distance = sc.sqrt((xs*xs) + (ys*ys))
	print centers[i], xc, yc, xs, ys, distance, i
	#lc, bc, rc = sc.mean(stuff[:,0]), sc.mean(stuff[:,1]), sc.mean(stuff[:,2])
	#ls, bs, rs = sc.std(stuff[:,0]), sc.std(stuff[:,1]), sc.std(stuff[:,2])
	#print lc, bc, rc
	#print ls, bs, rs
	#if i==36:
	#	for j in range(len(stuff[:,0])):		
	#		print X[j], Y[j], distance[j]
