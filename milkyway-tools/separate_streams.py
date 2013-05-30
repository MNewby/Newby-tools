import sys
sys.path.insert(0, '../utilities')
import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func
import astro_coordinates as coor
import glob

'''python script for quickly separating separation data
    Matthew Newby, March 23, 2012'''

def separation(folder="/home/newbym2/Desktop/sansSgr/"):
    """ Separates a separation file into distinct star files """
    files = glob.glob(folder+"/separation-[0-2][0-9]-new.txt")
    l = len(folder)
    for file in files:
        wedge = int(file[(l+11):(l+13)])
        if len(file) > 29:  name = file[(l+23):(l+29)]
        else:  name = "out"
        bg, s1, s2, s3 = [], [], [], []
        data = open(file, "r")
        for line in data:
            if int(line[0])==1:  s1.append(line[2:])
            elif int(line[0])==2:  s2.append(line[2:])
            elif int(line[0])==3:  s3.append(line[2:])
            else:  bg.append(line[2:])
        data.close()
        seps, tags, count = [bg, s1, s2, s3], ['bg-', 's1-', 's2-', 's3-'], 0
        for sep in seps:
            if sep != []:
                fileout = tags[count]+str(wedge)+name+"-new.txt"
                writefile = open(fileout, 'w')
                for line in sep:  writefile.write(line)
            count = count + 1
            writefile.close()
        print "{0}, {1} Separated".format(wedge, name)

def sep_lbr(folder1="/home/newbym2/Desktop/sansSgr/stars",
            folder2="/home/newbym2/Desktop/sansSgr/outfiles"):
    """ Matches an lbr starfile with an xyz starfile and saves the result"""
    files1 = glob.glob(folder1+"/stars-[0-2][0-9]*") #Star files
    files2 = glob.glob(folder2+"/out-*new*") #Separation files
    #files2 = glob.glob(folder2+"/*[7-8][2-9]*")
    #print files1
    #print files2
    dict1, dict2 = {}, {}
    for file in files1:
        temp_dict = {file : int(file[len(folder1)+7:len(folder1)+9])}
        #print temp_dict
        dict1.update(temp_dict)
    for file in files2:
        temp_dict = {file : int(file[len(folder2)+13:len(folder2)+15])}
        dict2.update(temp_dict)
    for file1 in files1:
        wedge = dict1[file1]
        # Test for matching file
        for item in dict2.items():
            if item[1]==wedge:  file2=item[0]; break
            else:  file2=None
        if file2==None:  print "# - No match for {0}".format(file1); continue
        print "# - Matching {0} and {1}".format(file1, file2)
        # Build new file
        flags, holder, count, skip = [], [], 0, 1
        data2 = open(file2, "r")
        for line in data2:
            flags.append(line[:2])
        data2.close()
        data1 = open(file1, "r")
        for line in data1:
            if skip==1: skip=0; continue  #skip the num_stars line
            holder.append(flags[count]+line)
            count = count + 1
        data1.close()
        fileout = "separation-"+str(wedge)+"-new.txt"
        writefile = open(fileout, "w")
        for line in holder:  writefile.write(line)
        writefile.close()

if __name__ == "__main__":
    #sep_lbr()
    separation()
    

print "# - All Done"
