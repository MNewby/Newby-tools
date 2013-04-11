import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt

'''python script for plotting gc data
Matthew Newby, August 1, 2010'''

x = 0
y = 2
x_label = 'Distance'
y_label = '$\mu$'
plot_file = 'results_HarrisCon_MUvDist.ps'
con_file = 'Convolve_results_Harris_cut.txt'
ori_file = 'original_results_table.txt'
label_file = 'label_file.txt'
con_list = []
if (con_file != ''):
    readfile = open(con_file, "r")
    for line in readfile:
        if (line[0] == "#"): continue
        if (line[0] == ""): continue
        con_list.append(list(line.split()))
    print '-Convolved data sucessfully loaded'
con_array = sc.array(con_list, dtype=float)
#print con_array
ori_list = []
if (ori_file != ''):
    readfile = open(ori_file, "r")
    for line in readfile:
        if (line[0] == "#"): continue
        if (line[0] == ""): continue
        ori_list.append(list(line.split()))
    print '-Original data sucessfully loaded'
ori_array = sc.array(ori_list, dtype=float)
#print ori_array
names = 0
label_list = []
if (label_file != ''):
    readfile = open(label_file, "r")
    for line in readfile:
        if (line[0] == "#"): continue
        if (line[0] == ""): continue
        if (names == 0):
            name_list = line.split(',')
            names = names + 1
        else: label_list.append(list(line.split()))
label_array = sc.array(label_list, dtype=float)
#print label_array, name_list
plt.figure()
#plt.plot(con_array[:,x], con_array[:,y], 'g-')
plt.scatter(con_array[:,x], con_array[:,y], c='blue', marker='+')
plt.scatter(ori_array[:,x], ori_array[:,y], c='green', marker='^')
for i in range(len(name_list)):
    plt.text(label_array[i][x], label_array[i][y], name_list[i])
plt.title('Plot of convolved data')
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.savefig(plot_file, papertype='letter')
plt.close('all')
print '-plot successfully created'