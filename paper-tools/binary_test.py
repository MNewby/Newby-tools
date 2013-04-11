import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func

'''python script for roughly studying maximum effect of binaries on CMD turnoffs.
Matthew Newby, Oct ??, 2011'''

import matplotlib.pyplot as plt

x1 = np.random.normal(0.25, 0.025, 1000)
x2 = sc.zeros(len(x1))

for i in range(len(x1)):
    if np.random.uniform() < 0.2:
        x2[i] = x1[i] + 0.05
    else:
        x2[i] = x1[i]

y1, z1 = np.histogram(x1, 40, range=(0.0, 0.4))
y2, z2 = np.histogram(x2, 40, range=(0.0, 0.4))

m1 = sc.zeros((len(y1), 2))
m1[:,0] = (z1[:-1]+0.005)
m1[:,1] = y1
m2 = sc.zeros((len(y2), 2))
m2[:,0] = (z2[:-1]+0.005)
m2[:,1] = y2
fi.write_data(m1, "ori_color.txt")
fi.write_data(m2, "new_color.txt")

plt.figure(1)
plt.bar(z2[:-1], y2, 0.01, color="blue")
plt.bar(z1[:-1],y1, 0.01, color="green", fill=None)
plt.show()



''' For magnitude study
""" For double-sided Gaussian """
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