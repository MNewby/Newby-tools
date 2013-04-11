import files as f
import scipy as sc

x = f.read_csv("6205_U_errors.csv")
l, w = x.shape
out = 0
for i in range(len(x[:,0])):
    if x[i,4] > 24.0:  out = out + 1
y = sc.zeros(((l-out),(w)), float)
j = 0
for i in range(len(x[:,0])):
    if x[i,4] <= 24.0:
        y[j,:] = x[i,:]
        j = j + 1
if j != (l-out): print 'CRAP!'
else: print (l-out), 'Stars in new set', l, 'stars in old'
f.write_csv(y, "6205_Uerrors_cut.csv")