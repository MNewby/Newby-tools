import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

'''python script for plotting u,g,r error fits on the same pane
Matthew Newby, May 4, 2011'''

#Get data
data1 = f.read_data('Pal5_Uerrors_cut2.csv', ',')
data2 = f.read_data('err_pal5_stars.txt')

#Assign data
u_data = data1[:,2]
u_err = data1[:,3]
g_data = data2[:,2]
g_err = data2[:,4]
r_data = data2[:,3]
r_err = data2[:,5]

#Initialize error functions
"""  Old Values
a_u = 0.02076735
b_u = 0.81309147
c_u = -19.61533299
a_g = 0.0
b_g = 0.790391
c_g = -19.8928
a_r = 0.0
b_r = 0.766309
c_r = -19.0334
"""
""" Best custom fit """
a_u = 2.71190687e-03   
b_u = 8.01759544e-01  
c_u = -1.92422351e+01
a_g = 3.10766606e-04   
b_g = 0.793349993
c_g = -19.9656922
a_r = -2.62791514e-05   
b_r = 0.797633236
c_r = -19.7509790

""" GNUplot fits 
a_u = 0.0    
b_u = 0.782646
c_u = -18.7956
a_g = 0.0
b_g = 0.791731
c_g = -19.9228
a_r = 0.0
b_r = 0.789554
c_r = -19.5491
"""

def error_function(x, params):
    y = params[0] + np.exp((params[1]*x)+params[2])
    return y

x_data = sc.arange(12.0, 24.0, 0.1)
y_u = error_function(x_data, [a_u, b_u, c_u])
y_g = error_function(x_data, [a_g, b_g, c_g])
y_r = error_function(x_data, [a_r, b_r, c_r])

# set data offsets
u_off = 1.0
g_off = 0.5
r_off = 0.0

plt.figure()
ax = plt.subplot(111)
#plot u data and fit
plt.errorbar(u_data, (u_err+u_off), fmt='o', ms=1, mfc='b', mec='b')  #data
plt.plot(x_data, (y_u+u_off), 'k-')  #Line fit
plt.plot([12.0, 24.0],[u_off, u_off],'k:')  # zero line
#plot g data and fit
plt.errorbar(g_data, (g_err+g_off), fmt='o', ms=1, mfc='g', mec='g')  #data
plt.plot(x_data, (y_g+g_off), 'k-')  #Line fit
plt.plot([12.0, 24.0],[g_off, g_off],'k:')  # zero line
#plot r data and fit
plt.errorbar(r_data, (r_err+r_off), fmt='o', ms=1, mfc='r', mec='r')  #data
plt.plot(x_data, (y_r+r_off), 'k-')  #Line fit
plt.plot([12.0, 24.0],[r_off, r_off],'k:')  # zero line

plt.xlim(12.0, 24.0)
plt.ylim(-0.1, 2.0)
ax.xaxis.set_minor_locator(MultipleLocator(0.4))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))

plt.xlabel('Apparent Magnitude', fontsize=10)
plt.ylabel('Magnitude Error', fontsize=10)
plt.text(12.5, 1.05, r'$\epsilon(u)+1.0$', fontsize=16)
plt.text(12.5, 0.55, r'$\epsilon(g)+0.5$', fontsize=16)
plt.text(12.5, 0.05, r'$\epsilon(r)$', fontsize=16)
plt.show()
print "# - Done"