import csv
import scipy as sc
import math as m
import numpy as np
import sys

'''This program reads files and searches for the 'results' tag:  $
2 argv:  in file and outfile
October 27, 2010'''

filename = sys.argv[1]
#filename = 'new_6205_runs_point10.txt'
if ( len(sys.argv) > 2 ):
    outfile = sys.argv[2]
else:
    outfile = 'new_results_out.txt'

readfile = open(filename, "r")
results = []
for line in readfile:
    if (line.strip() == ''): continue
    if (line[0] == "$"): results.append(line)
    elif (line[0] == '!'):  print line
    continue
out = []
for i in range(len(results)):
    line = results[i][1:].strip()
    out.append( line + '\n' )
writefile = open(outfile, 'w')
for i in range(len(out)):
    writestr = out[i]
    writefile.write(writestr)
writefile.close()
print '#-Results file written as', outfile