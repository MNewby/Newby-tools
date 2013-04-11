import files as f
import math as m
import numpy as np
import scipy as sc
import matplotlib
import matplotlib.pyplot as plt


rolloff_data = f.read_data("results_verif_cluscorr.txt")
sigmoid_data = f.read_data("results_verif_backcorr.txt")
straight_data = f.read_data("results_verif_nocorr.txt")
data_list = [rolloff_data, sigmoid_data, straight_data]

sigr = 2
errors = 6

out = []
for i in range(len(data_list)):
    data = data_list[i]
    series = []
    for j in range(5, 45, 2):
        holder = []
        weights = []
        for k in range(len(data[:,0])):
            if data[k,0] < j: continue
            elif data[k,0] >= (j+2): continue
            else:
                holder.append(data[k,sigr])
                weights.append( 1.0 / (data[k,errors]*data[k,errors]) )
        if (len(holder) == 0 ):  continue
        else:  series.append([(j+1), sc.average(holder, weights=weights )])
    out.append(sc.array(series))
#print out
fig = plt.figure()
ax1 = fig.add_subplot(111)
for i in range(len(data_list)):
    ax1.plot(out[i][:,0], out[i][:,1], 'k-')
x_tick = ax1.get_xticks()
ax1.set_xticks(x_tick)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
new_ticks = []
for i in range(len(x_tick)):
    new_ticks.append(round( (5.0*m.log10(x_tick[i]/0.01) + 4.18), 1) )
x_lim = ax1.get_xlim()
new_lim = (5.0*m.log10(x_lim[0]/0.01) + 4.18), (5.0*m.log10(x_lim[1]/0.01) + 4.18)
ax2 = ax1.twiny()
ax2.set_xlim(new_lim)
#ax2.set_xticks(new_ticks)
plt.xticks(x_tick, new_ticks)
plt.xticks(fontsize=10)
ax1.set_xlabel(r'$d_{eff}$')
ax2.set_xlabel(r'$\bar{g}_0$')
ax1.set_ylabel(r'$\sigma_r$')
plt.text(26.0, 1.32, "No detection loss", fontsize=10)
plt.text(26.0, 1.08, "With sigmoidal detection efficiency", fontsize=10, rotation=-9.0)
plt.text(26.0, 0.76, "With cluster detection efficiency", fontsize=10, rotation=-23.0)
plt.show()
print "#  Done"    