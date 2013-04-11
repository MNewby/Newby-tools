#! /usr/bin/env python

import csv
import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


"""This program plots data from the gc project

Matthew Newby (RPI), March 1, 2010
"""

#data, l,b

name = ['NGC 4147', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 5466', 'NGC 5904', 'NGC 6205', 'NGC 7078', 'NGC 7089', 'Pal 5', 'Sag Core']
l = [252.85, 332.96, 335.69, 42.21, 42.15, 3.86, 59.01, 65.01, 53.38, 0.85, 5.6]
b = [77.19, 79.76, 78.94, 78.71, 73.59, 46.80, 40.91, -27.31, -35.78, 45.86, -14.2]

plt.figure(1)
for i in range(len(l)):
    plt.scatter(l[i], b[i])
    plt.text(l[i], b[i], name[i])
plt.title('l and b plot of globular clusters')
plt.xlabel('l')
plt.ylabel('b')
plt.show()
