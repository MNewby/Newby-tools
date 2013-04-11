import scipy as sc
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import files as f

#Built for creating plots of isocrones and color-mag diagrams (in gnuplot)

d = raw_input('GC distance:')
distance = eval(d)
#isofile = raw_input('Isocrone file name:')
starfile = raw_input('Star file name:')

#x = f.read_data(isofile)
#l, w = x.shape
#a = sc.zeros(l)
#a = x[:,8] 
#b = sc.zeros(l)
#b = x[:,9]
#c = sc.zeros((l,2))
#c[:,0] = a
#c[:,1] = a - b
#if f.write_data(c, 'isocrone_xxxx_g-r.dat') == 1:
#    print '-Isocrone written succesfully'

y = f.read_data(starfile)
l, w = y.shape
j = sc.zeros(l)
j = y[:,2]
k = sc.zeros(l)
k = y[:,3]
m = sc.zeros((l,2))
m[:,1] = j - k
for i in range(l):
    j[i] = j[i] - 5.*(np.log10(distance*1000) - 1.)
m[:,0] = j

if f.write_data(m, 'xxxx_g-r.dat') == 1:
    print '-Isocrone written succesfully'


