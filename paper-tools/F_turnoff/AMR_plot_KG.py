import math as ma
from math import sin, cos
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func
import poly_fit as ply
#import matplotlib.gridspec as gridspec

def read_DA2005():
    data = []
    infile = open("AMR_DeAngeli2005.txt", 'r')
    for line in infile:
        if line.strip() == "":  continue
        if line.strip()[0] == "#":  continue
        holder = line.split()
        holder2 = [holder[-13], holder[-12], holder[-11], holder[-4], holder[-2],
                   holder[-3], holder[-1]]
        # Metal(ZW), Metal(CG), r(GC), Age(ZW), Age(CG), Age Err(ZW), Age Err(CG)
        for j in range(len(holder2)):
            holder2[j] = eval(holder2[j])
            if j==3 or j==5:  holder2[j] = 11.2*holder2[j]
            if j==4 or j==6:  holder2[j] = 10.9*holder2[j]
        data.append(holder2)
    infile.close()
    return sc.array(data)
    
def read_Dotter():
    data = []
    infile = open("AMR_Dotter2010.txt", 'r')
    for line in infile:
        if line.strip() == "":  continue
        if line.strip()[0] == "#":  continue
        holder = line.split()
        holder2 = [holder[-6], holder[-2], holder[-1], holder[-4]]
        # Metal (ZW), Age, age error, dist_mod
        for j in range(len(holder2)):
            holder2[j] = eval(holder2[j])
        if j==0:
            holder2[j] = 0.61 + 3.04*holder2[j] + 1.981*holder2[j]*holder2[j] + \
                0.532*holder2[j]*holder2[j]*holder2[j]
        data.append(holder2)
    infile.close()
    return sc.array(data)
    
def read_MF2009():
    data = []
    infile = open("AMR_MF2009.txt", 'r')
    for line in infile:
        if line.strip() == "":  continue
        if line.strip()[0] == "#":  continue
        holder = line.split()
        holder2 = [holder[-15], holder[-14], holder[-11], holder[-9],
                   holder[-10], holder[-8], holder[-1]]
        # Metal (ZW,CG), Age (ZW,CG), errors
        for j in range(len(holder2)):
            if j == 6:
                if holder2[j].strip()== "Old":  holder2[j] = 2
                elif holder2[j].strip()== "Young":  holder2[j] = 1
                else:  holder2[j] = 0
                continue
            holder2[j] = eval(holder2[j])
            if j>=2 and j<=5:  holder2[j] = 12.8*holder2[j]
        data.append(holder2)
    infile.close()
    return sc.array(data)
    
def read_Newby():
    data = []
    infile = open("AMR_Newby2011.txt", 'r')
    for line in infile:
        if line.strip() == "":  continue
        if line.strip()[0] == "#":  continue
        holder = line.split()
        holder2 = [holder[-9], holder[-8], holder[-7], holder[-3], holder[-1]]
        # Metal (CG), distance, Age, fit offset, error
        for j in range(len(holder2)):
            holder2[j] = eval(holder2[j])
        data.append(holder2)
    infile.close()
    return sc.array(data)

def read_clusters():
    data, holder = [], []
    infile = open("AMR_cluster_data.txt")
    for line in infile:
        if line.strip()[0] == '#':  continue
        if line.strip()[0] == '':  continue
        if line.strip()[0] == '*':
            if len(holder) > 0:
                data.append(holder)
                holder = []
            continue
        temp = line.split()
        for i in range(len(temp)):
            temp[i] = eval(temp[i].strip())
        holder.append(temp)
    if len(holder) > 0:
        data.append(holder)
    infile.close()
    return data

def clus_info(start, stop):
    # Data
    Newby = read_Newby()
    Dotter = read_Dotter()
    MF = read_MF2009()
    DA = read_DA2005()
    # cluster = [newby, dotter, MF, DA] -1=no data
    NGC4147 = [0, 8, 8, 8]
    NGC5024 = [1, 11, 11, 12]
    NGC5053 = [2, 12, 12, 50]
    NGC5272 = [3, 13, 14, 51]
    NGC5466 = [4, 15, 16, 52]
    NGC5904 = [5, 16, 17, 15]
    NGC6205 = [6, 24, 25, 58]
    NGC6341 = [7, 28, 29, 61]
    NGC7078 = [8, 49, 54, 38]
    NGC7089 = [9, 50, 55, 39]
    Pal5    = [10, -1, -1, -1]    
    index = sc.array([NGC4147, NGC5024, NGC5053, NGC5272, NGC5466, NGC5904,
                      NGC6205, NGC6341, NGC7078, NGC7089, Pal5])
    # Plotting - MAYBE PUT HARRIS, CG, OTHERS, ON HERE?
    shape = ['o', 's', '^', 'd', 'p', '>', '<', 'v', 'p', 'h']
    fc = ['black', 'blue', 'cyan', 'green', 'yellow', 'red', 'magenta', 'white', 'blue', 'red']
    ec = ['black', 'blue', 'cyan', 'green', 'yellow', 'red', 'magenta', 'black', 'black', 'black']
    fig = plt.figure(1)
    #gs = gridspec.Gridspec(4,1)
    plt.subplots_adjust(hspace=0.001)    
    """ Subplot 1 """
    #sp1 = plt.subplot2grid((4,1), (0,0), rowspan=2)
    #sp1 = plt.subplot(gs[0:1,0])
    for i in range(start,stop):  #skipping Pal 5
    ##    #separate by color (author), shape (cluster) <--switch!
        p = index[i]
        print "# - "
        # Newby
    #    plt.errorbar(Newby[p[0],2], Newby[p[0],0], xerr=(Newby[p[0],4]*2.0), marker='o', mfc=fc[i], mec=ec[i], ecolor=ec[i])
    #    plt.text(Newby[p[0],2], Newby[p[0],0], 'Newby', fontsize=6)
        print Newby[p[0],:]
    #    # Dotter
    #    plt.errorbar(Dotter[p[1],1], Dotter[p[1],0], xerr=Dotter[p[1],2], marker='^', mfc=fc[i], mec=ec[i], ecolor=ec[i])
    #    plt.text(Dotter[p[1],1], Dotter[p[1],0], 'Dotter', fontsize=6)
        print Dotter[p[1],:]
    #    # Marin-Franch
    #    plt.errorbar(MF[p[2],3], MF[p[2],1], xerr=MF[p[2],5], marker='s', mfc=fc[i], mec=ec[i], ecolor=ec[i])
    #    plt.text(MF[p[2],3], MF[p[2],1], 'Marin-Franch', fontsize=6)
        print MF[p[2],:]
    #    # DeAngeli
    #    plt.errorbar(DA[p[3],4], DA[p[3],1], xerr=DA[p[3],6], marker='d', mfc=fc[i], mec=ec[i], ecolor=ec[i])
    #    plt.text(DA[p[3],4], DA[p[3],1], 'DeAngeli', fontsize=6)
        print DA[p[3],:]
    #plt.ylim(-2.5,-0.5)
    #plt.ylabel(r"[Fe/H]")
    #plt.setp(sp1.get_xticklabels(), visible=False)
    """ Subplot 2 - Delta [Fe/H] vs Age"""
    #sp2 = plt.subplot2grid((4,1), (1,0))
    sp2 = plt.subplot(211)
    for i in range(start,stop):  #skipping Pal 5
        p = index[i]
        avg = sc.mean([Dotter[p[1],0], MF[p[2],1], DA[p[3],1]])
        std = sc.std([Dotter[p[1],0], MF[p[2],1], DA[p[3],1]])
        plt.errorbar((Newby[p[0],2]+Newby[p[0],4]), (avg-Newby[p[0],0]), marker='o')
    plt.plot([7.0,15.0], [0.0,0.0], 'k-')
    plt.plot([7.0,15.0], [std, std], 'k:')
    plt.plot([7.0,15.0], [-std, -std], 'k:')
    plt.ylabel(r"$\Delta [Fe/H]$")
    plt.setp(sp2.get_xticklabels(), visible=False)
    """ Subplot 3 - Delta Age vs Age"""
    sp3 = plt.subplot(212, sharex=sp2)
    #sp3 = plt.sp2 = plt.subplot2grid((4,1), (1,0))
    #plt.subplots_adjust(top=0.25, bottom=0.0)
    for i in range(start,stop):  #skipping Pal 5
        p = index[i]
        avg = sc.mean([Dotter[p[1],1], MF[p[2],3], DA[p[3],4]])
        std = sc.std([Dotter[p[1],1], MF[p[2],3], DA[p[3],4]])
        plt.errorbar((Newby[p[0],2]+Newby[p[0],4]), (avg-Newby[p[0],2]+Newby[p[0],4]),
           yerr=Newby[p[0],4], marker='o')
    plt.plot([7.0,15.0], [0.0,0.0], 'k-')
    plt.plot([7.0,15.0], [std, std], 'k:')
    plt.plot([7.0,15.0], [-std, -std], 'k:')
    plt.ylabel(r"$\Delta Age$")
    plt.xlabel(r"Age (Gyr)")
    #plt.xlim(7.0,15.0)
    plt.show()

def AMR_plot():
    # Data
    Newby = read_Newby()
    Dotter = read_Dotter()
    MF = read_MF2009()
    DA = read_DA2005()
    # Plotting
    plt.figure(1)
    #ax = plt.subplot(111)
    # AMR Trendline
    plt.plot([13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 6.0, 4.0, 1.5],
            [-1.7, -1.5, -1.35, -1.2, -1.05, -0.90, -0.65, -0.625, -0.7], 'k-', label="Muratov & Gnedin 2010")
    # Muratov and Gnedin blue box
    plt.fill([1.5, 13.0, 13.0, 7.5, 1.5],[-1.05, -1.05, -2.3, -1.38, -1.05], 'b', alpha=0.25)#, fill=None, hatch="/", color='blue')#, ls='dotted')
    #plt.plot([1.5, 13.0, 13.0, 7.5, 1.5],[-1.05, -1.05, -2.3, -1.38, -1.05], 'b-')
    # Muratov and Gnedin red box
    plt.fill([1.5, 1.5, 5.0, 6.8, 13.0],[-1.05, -0.35, 0.02, 0.02, -1.05], 'r', alpha=0.25)#, fill=None, hatch='\\', color='red') #ls='dotted')
    #plt.plot([1.5, 1.5, 5.0, 6.8, 13.0],[-1.05, -0.35, 0.05, 0.05, -1.05], 'r-')
    plt.errorbar(Dotter[:,1], Dotter[:,0], xerr=Dotter[:,2], fmt='o',
                 mec='black', mfc='black', ecolor='black', mew=1, markersize=5, elinewidth=0.5, label="Dotter et al. 2010")
    plt.errorbar(MF[:,3], MF[:,1], xerr=MF[:,5], fmt='s',
                 mec='black',mfc='black', ecolor='black', mew=1, markersize=5, elinewidth=0.5, label="Marin-Franch et al. 2009")
    plt.errorbar(DA[:,4], DA[:,1], xerr=DA[:,6], fmt='^',
                 mec='black',mfc='black', ecolor='black', mew=1, markersize=5, elinewidth=0.5, label="De Angeli at al. 2005")
    x,y,e = [], [], []
    for i in [0,1,3,4,5,6,7,9,10]:
        x.append(Newby[i,2]+Newby[i,3])
        y.append(Newby[i,0])
        e.append(Newby[i,4])
    plt.errorbar(x, y, xerr=e, fmt='o', mfc='cyan', mec='black', markersize=10, ecolor='black',
        elinewidth=1, label="Newby et al. 2011")
    x,y,e = [], [], []
    for i in [2,8]:
        x.append(Newby[i,2]+Newby[i,3])
        y.append(Newby[i,0])
        e.append(Newby[i,4])
    plt.errorbar(x, y, xerr=e, fmt='o', mfc='white', mec='cyan', mew=2.5, markersize=8, ecolor='black',
        elinewidth=1, label=None)
    plt.errorbar(x, y, xerr=None, fmt='o', mfc='white', alpha=0.5, mec='black', mew=1, markersize=10, label=None)
    plt.errorbar([12.0, 6.5],[-1.6, -0.65], xerr=None, fmt='o', mfc='yellow', mec='black', 
        mew=1.0, markersize=10, ecolor='black', elinewidth=1, label="This Paper")
    plt.xlabel(r"Age (Gyr)")
    plt.ylabel(r"[Fe/H]")
    plt.xlim(1.0, 15.0)
    plt.ylim(-2.5, 0.25)
    leg = plt.legend(loc='lower left', numpoints=1)
    for t in leg.get_texts():
        t.set_fontsize(8)
    plt.show()

def delta_plots():
    stuff = read_clusters()  # Stuff[cluster][paper][value]
    fig = plt.figure(1)
    plt.subplots_adjust(hspace=0.001)
    """ Get Newby Ages """
    base_a, base_e = [], []
    for i in range(10):  #All but Pal 5
        base_a.append(stuff[i][0][2])
        base_e.append(stuff[i][0][4])
    base_a.append(stuff[10][0][2])
    """ Distance plot - plot versus difference from Harris, DeAngeli/R-B """ #NEED ERROR BARS
    sp1 = plt.subplot(211)
    base_d, DA_d, Harris_d, dotter_d = [], [], [], []
    for i in range(10):   #skipping Pal 5
        base_d.append(stuff[i][0][1])
        DA_d.append(stuff[i][3][2] - base_d[i])
        Harris_d.append(stuff[i][4][0] - base_d[i])
        #dotter_hold = (10**((stuff[i][1][3]+1.0)/5.0) / 1000.0)
        dotter_hold = ( ( 10.**( (stuff[i][1][3])/5. ) )/ 100. )
        dotter_d.append(dotter_hold - base_d[i])
    #Append Pal 5 !!!
    #Adding Trendlines
    DA_fit = ply.poly_fit(0, base_a[:-1], DA_d, 0.1*sc.ones(len(DA_d)), 1, 0)
    DA_x = sc.arange(9.0, 15.0)
    DA_y = DA_x*DA_fit[1] + DA_fit[0]
    plt.plot(DA_x, DA_y, 'g:')
    Harris_fit = ply.poly_fit(0, base_a[:-1], Harris_d, 0.1*sc.ones(len(Harris_d)), 1, 0)
    Harris_x = sc.arange(9.0, 15.0)
    Harris_y = Harris_x*Harris_fit[1] + Harris_fit[0]
    plt.plot(Harris_x, Harris_y, 'k:')
    dotter_fit = ply.poly_fit(0, base_a[:-1], dotter_d, 0.1*sc.ones(len(dotter_d)), 1, 0)
    dotter_x = sc.arange(9.0, 15.0)
    dotter_y = dotter_x*dotter_fit[1] + dotter_fit[0]
    plt.plot(dotter_x, dotter_y, 'b:')
    plt.plot([9.0,14.0], [0.0,0.0], 'k--')
    plt.errorbar(base_a[:-1], DA_d, yerr=None, fmt='^', mec='green', mfc='green', label="De Angeli et al. 2005")
    plt.errorbar(base_a[:-1], Harris_d, yerr=None, fmt='*', mec='black', mfc='black', label="Harris Catalog (2010)")
    plt.errorbar(base_a[:-1], dotter_d, yerr=None, fmt='o', mec='blue', mfc='blue', label="Dotter et al. 2010")
    plt.setp(sp1.get_xticklabels(), visible=False)
    plt.ylabel(r'$\Delta$ Distance (kpc)')
    plt.ylim(-5.0, 5.0)
    leg = plt.legend(loc='upper left')
    for t in leg.get_texts():
        t.set_fontsize(8)
    """ Age plot - plot versus difference from Harris, DeAngeli/R-B """
    sp2 = plt.subplot(212, sharex=sp1)
    dotter_a, MF_a, DA_a = [], [], []
    dotter_e, MF_e, DA_e = [], [], []
    for i in range(10):  #No age for Pal 5
        dotter_a.append(stuff[i][1][1] - base_a[i])
        dotter_e.append(stuff[i][1][2])
        MF_a.append(stuff[i][2][3] - base_a[i])
        MF_e.append(stuff[i][2][5])
        DA_a.append(stuff[i][3][4] - base_a[i])
        DA_e.append(stuff[i][3][6])
    plt.plot([9.0,14.0], [0.0,0.0], 'k--')
    #Adding Trendlines
    dotter_fit = ply.poly_fit(0, base_a[:-1], dotter_a, dotter_e, 1, 0)
    dotter_x = sc.arange(9.0, 15.0)
    dotter_y = dotter_x*dotter_fit[1] + dotter_fit[0]
    plt.plot(dotter_x, dotter_y, 'b:')
    MF_fit = ply.poly_fit(0, base_a[:-1], MF_a, MF_e, 1, 0)
    MF_x = sc.arange(9.0, 15.0)
    MF_y = MF_x*MF_fit[1] + MF_fit[0]
    plt.plot(MF_x, MF_y, 'r:')
    DA_fit = ply.poly_fit(0, base_a[:-1], DA_a, DA_e, 1, 0)
    DA_x = sc.arange(9.0, 15.0)
    DA_y = DA_x*DA_fit[1] + DA_fit[0]
    plt.plot(DA_x, DA_y, 'g:')
    #plt.errorbar((sc.array(base_a[:-1]) - 0.05), sc.zeros(len(base_a[:-1])), yerr=base_e, fmt=None, ecolor='black')
    plt.errorbar(base_a[:-1], dotter_a, yerr=dotter_e, fmt='o', mec='blue', mfc='blue',
                 ecolor='blue', label="Dotter et al. 2010")
    plt.errorbar(base_a[:-1], MF_a, yerr=MF_e, fmt='s', mec='red', mfc='red',
                 ecolor='red', label="Marin-Franch et al. 2009")
    plt.errorbar(base_a[:-1], DA_a, yerr=DA_e, fmt='^', mec='green', mfc='green',
                 ecolor='green', label="De Angeli et al. 2005")
    plt.ylabel(r"$\Delta$ Age (Gyr)")
    plt.xlabel("Age (Gyr)")
    plt.xlim(9.0,14.0)
    plt.ylim(-5.0, 5.0)
    plt.yticks(sc.arange(-4.0,4.1,1.0))
    leg = plt.legend(loc='lower left', numpoints=1)
    for t in leg.get_texts():
        t.set_fontsize(8)
    plt.show()
    return 0
    
if __name__ == "__main__":
    #delta_plots()
    #read_clusters()[0]
    AMR_plot()
    #clus_info(0,10)
