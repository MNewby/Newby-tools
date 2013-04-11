import math as m
import numpy as np
import scipy as sc
import files as f
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt

'''python script for stacking histogram fits to convolved globular cluster data
Matthew Newby, May 19, 2011'''

# column 1 = M_g;  column 2= bin height
near = f.read_data('5272_10.txt')
far = f.read_data('5272_44.txt')
near_params = [4.3264132743867521, 0.39067339740578616, 0.76703139246245022, 159.36335049498547]
far_params = [3.9715809742388561, 0.34668679144488945, 0.55005316028356843, 38.20385043203013]
step_size = 0.2
half_step = step_size / 2.0
edge = step_size * 2.0

fig = plt.figure()
plt.subplots_adjust(hspace=0.001)

subs=[]
""" Near (top) plot """
sp = plt.subplot(211)
subs.append(sp)
#Make histogram
x_near, y_near = [], []
for i in range(len(near[:,0])):
    x_near.append(near[i,0]-half_step)
    x_near.append(near[i,0]+half_step)
    y_near.append(near[i,1]-half_step)
    y_near.append(near[i,1]+half_step)
#Make Gaussian fit line
x_line = sc.arange( (np.ma.min(near[:,0]) - edge), (np.ma.max(near[:,0]) + edge), 0.01)
y_line = sc.zeros(len(x_line))
for i in range(len(x_line)):
    if x_line[i] < near_params[0]:
        sigma = near_params[1]
    else:
        sigma = near_params[2]
    exponent = -0.5 * ( (x_line[i] - near_params[0]) / sigma )**2
    stuff = m.exp(exponent)
    y_line[i] = stuff*near_params[3]
plt.plot(x_near, y_near, 'k-')
plt.plot(x_line, y_line, 'g-')
plt.text(1.0, 120.0, "10.4 kpc")
plt.ylabel("counts", fontsize=10)
plt.ylim(0.0, 160.0)
plt.yticks(sc.arange(20.0, 180.0, 20.0))
plt.setp(sp.get_xticklabels(), visible=False)

""" Far (bottom) plot """
sp = plt.subplot(2,1,2, sharex=subs[0])
subs.append(sp)
#Make histogram
x_far, y_far = [], []
for i in range(len(far[:,0])):
    x_far.append(far[i,0]-half_step)
    x_far.append(far[i,0]+half_step)
    y_far.append(far[i,1]-half_step)
    y_far.append(far[i,1]+half_step)
#Make Gaussian fit line
x_line = sc.arange( (np.ma.min(far[:,0]) - edge), (np.ma.max(far[:,0]) + edge), 0.01)
y_line = sc.zeros(len(x_line))
for i in range(len(x_line)):
    if x_line[i] < far_params[0]:
        sigma = far_params[1]
    else:
        sigma = far_params[2]
    exponent = -0.5 * ( (x_line[i] - far_params[0]) / sigma )**2
    stuff = m.exp(exponent)
    y_line[i] = stuff*far_params[3]
plt.plot(x_far, y_far, 'k-')
plt.plot(x_line, y_line, 'g-')
plt.text(1.0, 40.0, "44.0 kpc")
plt.ylabel("counts")
plt.xlabel(r"$M_g$")
plt.ylim(0.0, 50.0)
plt.xlim(0.0, 9.0)


""" Finish up """
plt.show()
print "#---Done"