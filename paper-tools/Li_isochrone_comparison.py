import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func
import astro_coordinates as coor

'''python script for quickly analyzing data.
Matthew Newby, November 15,2010'''

data = fi.read_data("Grillmair_2011_isochrone.dat")
data2 = fi.read_data("isochrone_8gyr_n1p0.dat")
data3 = fi.read_data("isochrone_7gyr_n0p8.dat")
data4 = fi.read_data("isochrone_6gyr_n0p8.dat")

g = coor.getg(9.7, data[:,8])
g_minus_r = data[:,8] - data[:,9]
g_minus_i = data[:,8] - data[:,10]

g2 = coor.getg(9.7, data2[:,8]) #+0.43
g3 = coor.getg(9.7, data3[:,8]) #+0.63
g4 = coor.getg(9.7, data4[:,8]) #+0.63

s1 = 0.41
s2 = 0.61
s3 = 0.52

M1 = 15.0 #coor.getM(15.0, 9.7)
M2 = 19.5 #coor.getM(19.5, 9.7)
M3 = 22.0 #coor.getM(22.0, 9.7)

plt.figure(1)
plt.subplots_adjust(wspace=0.001)

sp1 = plt.subplot(121)
plt.plot(g_minus_r, g, "k-", label="Grillmair")
#plt.plot((data2[:,8]-data2[:,9]), g2, "r:", label="__nolabel__")
#plt.plot((data2[:,8]-data2[:,9]), (g2+s1), "r--", label="8.0 Gyr, [Fe/H]=-1.0")
plt.plot((data3[:,8]-data3[:,9]), g3, "b:", label="__nolabel__")
plt.plot((data3[:,8]-data3[:,9]), (g3+s2), "b--", label="7.0 Gyr, [Fe/H]=-0.8")
plt.plot((data4[:,8]-data4[:,9]), g4, "g:", label="__nolabel__")
plt.plot((data4[:,8]-data4[:,9]), (g4+s3), "g--", label="6.0 Gyr, [Fe/H]=-0.8")
#plt.plot([-0.5, 1.5], [M1,M1], "k:")
plt.plot([-0.5, 1.5], [M2,M2], "k:")
plt.plot([-0.5, 1.5], [M3,M3], "k:")
#plt.ylim(8.0, 2.0)
plt.xlim(-0.0, 0.8)
plt.xlabel(r"$g-r$")
plt.ylabel(r"$g$")
#subby 2
sp2 = plt.subplot(122, sharey=sp1)
plt.plot(g_minus_i, g, "k-", label="Grillmair 2011")
#plt.plot((data2[:,8]-data2[:,10]), g2, "r:", label="__nolabel__")
#plt.plot((data2[:,8]-data2[:,10]), (g2+s1), "r--", label="8.0 Gyr, [Fe/H]=-1.0")
plt.plot((data3[:,8]-data3[:,10]), g3, "b:", label="__nolabel__")
plt.plot((data3[:,8]-data3[:,10]), (g3+s2), "b--", label="7.0 Gyr, [Fe/H]=-0.8")
plt.plot((data4[:,8]-data4[:,10]), g4, "g:", label="__nolabel__")
plt.plot((data4[:,8]-data4[:,10]), (g4+s3), "g--", label="6.0 Gyr, [Fe/H]=-0.8")
plt.plot([-0.5, 1.5], [M2,M2], "k:")
plt.plot([-0.5, 1.5], [M3,M3], "k:")
plt.ylim(24.1, 16.9)
plt.xlim(-0.0, 0.8)
plt.xlabel(r"$g-i$")
leg = plt.legend(loc='lower left')
for t in leg.get_texts():
    t.set_fontsize(8)
plt.show()

"""
#RA:  30.6863 : 30.7806
#DEC:  -3.29063 : -3.19688
params = [4.11427748, 0.57501995, 0.8918551, 8.91381507]
exp_params = [4.18, 0.36, 0.89, 9.0]

data = fi.read_data("whiting1_clus.csv", ",")
#data = fi.read_data("whiting1ug_kgrabows.csv", ",") 
g0 = data[:,2]
Mg = coor.getM(g0, 30.1)  #Uses distance fomula and Harris distance to get absolute magnitudes

hist, edges = np.histogram(Mg, 40, range=(0.0, 8.0))
plt.figure(1)
plt.bar(edges[:-1], hist, 0.2, facecolor='blue') #fill=False)
x = sc.arange(0.0, 8.0, 0.1)
plt.plot(x, func.get_2gauss_y(x, params), 'r', lw=2)
plt.plot(x,func.get_2gauss_y(x, exp_params), 'r--', lw=1)
plt.xlabel(r"$M_g$")
plt.show()
"""

#data_out = sc.zeros((len(hist),2), float)
#for i in range(len(hist)):
#    data_out[i,0] = hist[i]
#    data_out[i,1] = edges[i] + 0.1
#fi.write_data(data_out, "Whiting1_hist.txt")

"""
# l,b,dered_u,dered_g,dered_r,dered_i,dered_z

data1 = fi.read_data("NGC_5024_clean.csv", ",")
data2 = fi.read_data("NGC_5024_semi_clean.csv", ",")
data3 = fi.read_data("NGC_5024_dirty.csv", ",")

plt.figure(1)
y1, x1 = data1[:,3], (data1[:,3]-data1[:,4])
plt.scatter(x1, y1, 1, 'k', 'o')
plt.ylim(26.0, 15.0)
plt.xlim(-0.3, 0.6)
plt.ylabel(r'$g_0$')
plt.xlabel(r"$(g-r)_0$")
plt.text(-0.2, 16.0, "NGC 5024 - clean", fontsize=8)

plt.figure(2)
y2, x2 = data2[:,3], (data2[:,3]-data2[:,4])
plt.scatter(x2, y2, 1, 'k', 'o')
plt.ylim(26.0, 15.0)
plt.xlim(-0.3, 0.6)
plt.ylabel(r'$g_0$')
plt.xlabel(r"$(g-r)_0$")
plt.text(-0.2, 16.0, "NGC 5024 - semi-clean", fontsize=8)

plt.figure(3)
y3, x3 = data3[:,3], (data3[:,3]-data3[:,4])
plt.scatter(x3, y3, 1, 'k', 'o')
plt.ylim(26.0, 15.0)
plt.xlim(-0.3, 0.6)
plt.ylabel(r'$g_0$')
plt.xlabel(r"$(g-r)_0$")
plt.text(-0.2, 16.0, "NGC 5024 - dirty", fontsize=8)

plt.show()
plt.close('all')
"""


""" for randomly selecting stars from the unmatched star file 
import astro_coordinates as coor

data = fi.read_data("no_matches_82.txt", ",")
samples = np.random.uniform(0, len(data[:,0]), 20)
for i in range(len(samples)):
    l, b = data[samples[i],0], data[samples[i],1]
    ra, dec = coor.lbToEq(l,b)
    print float(ra), float(dec)
"""

'''
""" For double-sided Gaussian """
import matplotlib.pyplot as plt
N1 = 642 #N2*(0.36/0.76)
N2 = 1357 #2000.0 / ((0.36/0.76)+1)

a1 = np.random.normal(4.2, 0.36, N1)
a2 = np.random.normal(4.2, 0.76, N2)
holder = []
for i in range(len(a1)):
    if a1[i] < 4.2:  holder.append(a1[i])
for i in range(len(a2)):
    if a2[i] > 4.2:  holder.append(a2[i])
x1 = sc.array(holder)
#END OF DOBLE_GAUSSIAN

#x1 = np.random.normal(4.2, 0.8, 1000)
x2 = sc.zeros(len(x1))
for i in range(len(x1)):
    if np.random.uniform() < 0.2:
        x2[i] = x1[i] - 0.75
    else:  x2[i] = x1[i]
y1, z1 = np.histogram(x1, 40, range=(0.0, 8.0))
y2, z2 = np.histogram(x2, 40, range=(0.0, 8.0))

m1 = sc.zeros((len(y1), 2))
m1[:,0] = (z1[:-1]+0.1)
m1[:,1] = y1
m2 = sc.zeros((len(y2), 2))
m2[:,0] = (z2[:-1]+0.1)
m2[:,1] = y2
fi.write_data(m1, "ori_dhist.txt")
fi.write_data(m2, "new_dhist.txt")

b1, c1 = np.histogram(a1, 40, range=(0.0, 8.0))
b2, c2 = np.histogram(a2, 40, range=(0.0, 8.0))

plt.figure(1)
plt.bar(z1[:-1],y1, 0.2, color="green")
plt.bar(c1[:-1], b1, 0.2, edgecolor="red", fill=None)
plt.bar(c2[:-1], b2, 0.2, edgecolor="blue", fill=None)
#plt.bar(z2[:-1], y2, 0.2, color="green")
plt.show()
'''
"""
data = open("AMR_Newby2011.txt", "r")
for line in data:
    if line[0]=="#":  continue
    holder = line.split()
    raw_age = eval(holder[7])
    offset = eval(holder[11])
    error = eval(holder[13])
    print "{0} $\pm$ {1}".format(round((raw_age+offset),2),round(error,2))
"""
"""
iso = fi.read_data("./Girardi_Isochrones/5904_0.dat")
data = fi.read_data("noU_NGC_5904_cluster.csv", ",")

plt.figure()
plt.scatter((data[:,4]-data[:,2]), data[:,4]-14.515, 1, 'k')
#plt.plot((iso[:,7]-iso[:,8]), (iso[:,7] - 5.*(sc.log10(8.0*1000) - 1.)), 'r-')
plt.plot((iso[:,7]-iso[:,8]), (iso[:,7]), 'r-')
plt.xlabel(r"$(u-g)_0$")
plt.ylabel(r"$u_0$")
plt.ylim(10.0, -5.0)
plt.show()
"""

"""
def make_ticklabels_invisible(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_visible(False)


plt.figure(0)
ax1 = plt.subplot2grid((3,3), (0,0), colspan=3)
ax2 = plt.subplot2grid((3,3), (1,0), colspan=2)
ax3 = plt.subplot2grid((3,3), (1, 2), rowspan=2)
ax4 = plt.subplot2grid((3,3), (2, 0))
ax5 = plt.subplot2grid((3,3), (2, 1))

plt.suptitle("subplot2grid")
make_ticklabels_invisible(plt.gcf())
plt.show()
"""

"""
d0 = fi.read_data("./Girardi_Isochrones/5024_0.dat")
d1 = fi.read_data("./Girardi_Isochrones/5024_1.dat")
d2 = fi.read_data("./Girardi_Isochrones/5024_2.dat")
d3 = fi.read_data("./Girardi_Isochrones/5024-1.dat")
d4 = fi.read_data("./Girardi_Isochrones/5024-2.dat")

plt.figure(1)
plt.plot(d0[:,8]-d0[:,9], d0[:,8])
plt.plot(d1[:,8]-d1[:,9], d1[:,8])
plt.plot(d2[:,8]-d2[:,9], d2[:,8])
plt.plot(d3[:,8]-d3[:,9], d3[:,8])
plt.plot(d4[:,8]-d4[:,9], d4[:,8])
plt.xlim(0.0,0.4)
plt.ylim(10.0, 0.0)
plt.show()

import rotation3D as r3d

print r3d.roty(-90.0)
print r3d.roty(0.0)
print r3d.rotz(0.0)

print r3d.roty(-90.0)*sc.matrix([[1.0, 1.0, 1.0]]).T
print (r3d.roty(-90.0)*r3d.roty(0.0)*r3d.rotz(0.0)).I*sc.matrix([[1.0, 1.0, 1.0]]).T
"""

"""
data = fi.read_data("2sig_nothresh.txt")
new_data = []
g = lambda x:  (8.0/2.5)*x + 21.2
h = lambda x:  (-22.0)*x + 26.0
for i in range(len(data[:,0])):
    if data[i,2] < 18.0:
        if data[i,2] < h(data[i,2]-data[i,3]):  continue
        if data[i,2] > h(data[i,2]-data[i,3]-0.2):  continue
    if data[i,2] > 22.0:
        if data[i,2] > g(data[i,2]-data[i,3]):  continue
    new_data.append(data[i,:])
new_data = sc.array(new_data)
g_minus_r = new_data[:,2] - new_data[:,3]
fig = plt.figure()
plt.scatter(g_minus_r, new_data[:,2], 1, 'k', 'o')
plt.ylim(25.0, 15.0)
plt.xlim(-1.0, 2.0)
plt.ylabel(r'$g_0$')
plt.xlabel(r"$(g-r)_0$")
plt.show()
fname = raw_input("save data? [filename to save, empty to skip]")
if fname != "":
    fi.write_data(new_data, fname, header="NGC_6205")
#(-1.0, 18.0; 1.5, 26.0)
#(0.0, 26.0; 0.5, 15.0)
"""