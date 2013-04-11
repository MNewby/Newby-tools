import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib as mpl
import matplotlib
import matplotlib.pyplot as plt

'''python script for SEGUE plate plotting
Matthew Newby, July 22, 2010'''

file1 = "SEGUE_plates.dat"
#file2 = "SEGUE_test.dat"
data1 = f.read_data(file1)
#data2 = f.read_data(file2)
l1, w1 = data1.shape
#l2, w2 = data2.shape()

in_disk = 0
plt.figure()
for i in range (l1):
    if (abs(data1[i,5]) <= 20.0):
        plt.scatter(data1[i,4], data1[i,5], c='blue', marker='o')
        print str(data1[i,0]), "is in the disk, b =", data1[i,5]
        in_disk = in_disk + 1
    else:
        plt.scatter(data1[i,4], data1[i,5], c='red', marker='x')
    if (data1[i,2] < -1000.0):  faint_id = "---"
    else:  faint_id = str(data1[i,2])
    plate_id = str(data1[i,0]) + "/" + faint_id
    plt.text(data1[i,4], data1[i,5], (plate_id), fontsize=8)
print in_disk, 'plates in the galactic disk, out of', l1, 'total plates'
plt.title('Main SEGUE plates')
plt.xlabel('l')
plt.ylabel('b')
plt.show()
print '-process complete'