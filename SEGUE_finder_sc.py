import math as m
import numpy as np
import scipy as sc
import files as f

'''python script for finding SEGUE plate ra, dec.
Matthew Newby, Aug 16, 2010'''

data = f.read_csv('SEGUE_t2_mnewby.csv')
#plate,fiberid,ra,dec,l,b,feha,teffa,logga,alphafe,elodierv
length, width = data.shape
#finding number and ids of unique plates
plate_id = []
for i in range(length):
    if (plate_id.count(data[i,0])==0):
        plate_id.append(data[i,0])
print len(plate_id), plate_id
#finding min and max ra and dec for each plate
plate_values = []
for i in range(len(plate_id)):
    field = [1000, -1000, 1000, -1000]
    for j in range(length):
        if (data[j,0]==plate_id[i]):
            if (data[j,2] < field[0]):  field[0]=int(data[j,2])
            if (data[j,2] > field[1]):  field[1]=(int(data[j,2])+1)
            if (data[j,3] < field[2]):  field[2]=int(data[j,3])
            if (data[j,3] > field[3]):  field[3]=(int(data[j,3])+1)
    plate_values.append(field)
    #print '(( ra BETWEEN', field[0], 'AND', field[1], \
    #') AND ( dec BETWEEN', field[2], 'AND', field[3], ')) OR'
print 'All plates searched'
plate_values.sort()
print plate_values
kill_index = []
for i in range((len(plate_values)-1)):
    for j in range(i+1, len(plate_values)):
        if (plate_values[i][0] == plate_values[j][0]):
            if (plate_values[i][1] == plate_values[j][1]):
                if (plate_values[i][2] == plate_values[j][2]):
                    if (plate_values[i][3] == plate_values[j][3]):
                        kill_index.append(j)
print 'Number of duplicate plates:', len(kill_index)
kill_index.sort()
for i in range(len(kill_index)):
    print kill_index[i], plate_values[kill_index[i]]
    plate_values.pop(kill_index[i])
    for j in range(len(kill_index)):  kill_index[j] = kill_index[j]-1
for i in range(len(plate_values)):
    print '(( ra BETWEEN', plate_values[i][0], 'AND', plate_values[i][1], \
    ') AND ( dec BETWEEN', plate_values[i][2], 'AND', plate_values[i][3], ')) OR'
print 'Done'