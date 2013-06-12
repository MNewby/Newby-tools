#import sys
#sys.path.insert(0, '../utilities')
#import functions as func
import math as ma
import numpy as np
import scipy as sc
import matplotlib
import matplotlib.pyplot as plt


file1 = "/home/newbym2/Dropbox/Research/Cetus_Polar/bhb_histogram.txt"
#file1 = "/home/newbym2/Dropbox/Research/Cetus_Polar/bhb_histogram_wide.txt"

data1 = sc.array(sc.loadtxt(file1, delimiter=","))

Xedges = data1[:,0] - 5.0
colX = sc.linspace(80,180,150)
xlocs = sc.arange(80.0, 180.1, 20.0)
xlabels = [r"$80$", r"$100$", r"$120$", r"$140$", r"$160$", r'$180$']
#for lab in xlocs:  xlabels.append(r"$"+str(lab)+"$")
ylocs = range(0, 21, 4)
ylabels = []
for lab in ylocs:  ylabels.append(r"$"+str(lab)+"$")
xname = '$l$ $(\hspace{-0.15}^\circ\hspace{-0.20})$'

xt, yt = 100.0, 18.0 # These set the position of the text in each subplot

bar1 = data1[:,1]
bar2 = data1[:,2]
bar3 = data1[:,3]
bar4 = data1[:,4]
bar5 = data1[:,5]

# This is a really ugly way to write out Gaussians...
Y1 = 3.4817747477437186**2*sc.exp(-(colX-140.7302048737979)**2/2/11.620129178363934**2)+0.78712169834133416
Y2 = 3.7244964876390734**2*sc.exp(-(colX-136.68217644170682)**2/2/5.5057461891887112**2)+2.6822071366004825
Y3 = 3.5877437862458397**2*sc.exp(-(colX-149.96711235113372)**2/2/15.107198276605352**2)+0.48085241051210986
Y4 = 3.9275520725880302**2*sc.exp(-(colX-142.19516947043576)**2/2/4.9711557834354672**2)+2.7128320752226713
Y5 = 3.4719517680362308**2*sc.exp(-(colX-144.74349790836359)**2/2/9.7924634139585454**2)+3.3554256983677386

#func.gauss_plus_floor()


fig = plt.figure(1)
plt.subplots_adjust(hspace=0.1, wspace=0.001)

sp4 = plt.subplot(234)
sp4.bar(Xedges, bar4, width=10.0, color="grey", hatch="\\\\")
sp4.bar(Xedges[:-3], bar4[:-3], width=10.0, color="grey", hatch=None)
sp4.plot(colX, Y4, "k-")
sp4.text(xt,yt, r'$-50^\circ < b < -42^\circ$')
plt.xlim(80.0, 180.0)
plt.xlabel(xname, fontsize = 14)
plt.xticks(xlocs, xlabels, fontsize=10)
plt.yticks(ylocs, ylabels, fontsize=10)
#plt.ylim(0.0, 20.0)
plt.ylabel(r"$N$", fontsize = 14)

sp1 = plt.subplot(231, sharex=sp4)
sp1.bar(Xedges, bar1, width=10.0, color="grey", hatch="\\\\")
sp1.bar(Xedges[:-3], bar1[:-3], width=10.0, color="grey", hatch=None)
sp1.plot(colX, Y1, "k-")
sp1.text(xt,yt, r'$-70^\circ < b < -62^\circ$')
#plt.setp(sp1.get_xticklabels(), visible=False)
#plt.xlim(80.0, 180.0)
#plt.xlabel('$l$ ($^\circ$)', fontsize = 14)
plt.xticks(xlocs, xlabels, fontsize=10)
plt.ylim(0.0, 20.0)
plt.ylabel(r"$N$", fontsize = 14)
plt.yticks(ylocs, ylabels, fontsize=10)

sp5 = plt.subplot(235, sharey=sp4)
sp5.bar(Xedges, bar5, width=10.0, color="grey", hatch="\\\\")
sp5.bar(Xedges[:-3], bar5[:-3], width=10.0, color="grey", hatch=None)
sp5.plot(colX, Y5, "k-")
sp5.text(xt,yt, r'$-42^\circ < b < -30^\circ$')
plt.setp(sp5.get_yticklabels(), visible=False)
plt.xlim(80.0, 180.0)
plt.xlabel(xname, fontsize=14)
plt.xticks(xlocs[1:], xlabels[1:], fontsize=10)
plt.ylim(0.0, 20.0)
#plt.ylabel(r"$N$", fontsize = 14)
#plt.yticks(ylocs, ylabels, fontsize=10)

sp2 = plt.subplot(232, sharex=sp5, sharey=sp1)
sp2.bar(Xedges, bar2, width=10.0, color="grey", hatch="\\\\")
sp2.bar(Xedges[:-3], bar2[:-3], width=10.0, color="grey", hatch=None)
sp2.plot(colX, Y2, "k-")
sp2.text(xt,yt, r'$-62^\circ < b < -56^\circ$')
#plt.setp(sp2.get_xticklabels(), visible=False)
plt.setp(sp2.get_yticklabels(), visible=False)
#plt.xlim(80.0, 180.0)
#plt.xlabel(xname, fontsize=14)
plt.xticks(xlocs[1:], xlabels[1:], fontsize=10)
#plt.ylim(0.0, 20.0)
#plt.ylabel(r"$N$", fontsize = 14)
#plt.yticks(ylocs, ylabels, fontsize=10)

sp3 = plt.subplot(233, sharey=sp1)
sp3.bar(Xedges, bar3, width=10.0, color="grey", hatch="\\\\")
sp3.bar(Xedges[:-3], bar3[:-3], width=10.0, color="grey", hatch=None)
sp3.plot(colX, Y3, "k-")
sp3.text(xt,yt, r'$-56^\circ < b < -50^\circ$')
plt.setp(sp3.get_yticklabels(), visible=False)
plt.xlabel(xname, fontsize = 14)
plt.xticks(xlocs[1:], xlabels[1:], fontsize=10)
plt.xlim(80.0, 180.0)
plt.ylim(0.0, 20.0)
#plt.ylabel(r"$N$", fontsize = 14)
#plt.yticks(ylocs, ylabels, fontsize=10)

plt.show()

