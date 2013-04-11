import math as ma
import numpy as np
import scipy as sc
import files as fi

arcsec = 0.000278

# modify!  Make test data sets!  Clean vs unclean vs coadd - contamination, completeness.

class Compare:
    def __init__(self, data1, data2, c1=0, c2=0, sub=100):
        """ c1 and c2 are columns to sort, sub is fineness of sort """
        # Sort data
        self.data1 = data1[data1[:,c1].argsort(),:]
        print "# Data 1 sorted"
        self.data2 = data2[data2[:,c2].argsort(),:]
        print "# Data 2 sorted"
        self.l1 = len(data1[:,0])
        self.l2 = len(data2[:,0])
        # Flags for data; 
        self.flags1 = sc.zeros(self.l1, int)
        self.flags2 = sc.zeros(self.l2, int)
        print "# Flags initialized"
        # Index data
        self.index1 = create_indices(self.data1, sub, c1)
        self.index2 = create_indices(self.data2, sub, c2)
        print "# Indices Generated"
    
def match_stars(sort, sep=arcsec):
    """ columns are l and b or ra and dec; matches according to sep
        flags:  0=unknown, 1=matched """
    for i in range(sort.l1):
        if sort.flags1[i] > 0: continue
        l0, b0 = sort.data1[i,0], sort.data1[i,1]
        k=0  #Find the first index that overshoots l0, search from the index prior to that
        try:
            while sort.index2[k][0] <= l0: k=k+1
            a,b = sort.index2[k-1][1], sort.index2[k][1]
        except IndexError:
            a,b = sort.index2[-1][1], (sort.l2-1)
        for j in range(a,b):
            #if sort.flags2[j] > 0: continue  #skips matching if matched, but misses double matches
            l1, b1 = sort.data2[j,0], sort.data2[j,1]
            d = ma.sqrt((l1-l0)*(l1-l0) + (b1-b0)*(b1-b0))
            if d < sep:
                sort.flags1[i], sort.flags2[j] = 1, 1
                # NO BREAK HERE - this way, multiple matches will be taken care of.
        if i % 10000 == 0:
            print "# - Searched {0} stars, {1} matches".format(i, sc.sum(sort.flags1))
    print "# --- Search of {0} stars complete, {1} matches".format(i, sc.sum(sort.flags1))

def separate_stars(sort, files=["matched.txt", "nomatch1.txt", "nomatch2.txt"]):
    # Matched stars
    holder = []
    for i in range(sort.l1):
        if sort.flags1[i] == 1:
            holder.append(sort.data1[i,:])
    fi.write_data(sc.array(holder), files[0])
    s1 = len(holder)
    # Unmatched 1
    holder = []
    for i in range(sort.l1):
        if sort.flags1[i] == 0:
            holder.append(sort.data1[i,:])
    fi.write_data(sc.array(holder), files[1])
    s2 = len(holder)
    # Unmatched 2
    holder = []
    for i in range(sort.l2):
        if sort.flags2[i] == 0:
            holder.append(sort.data2[i,:])
    fi.write_data(sc.array(holder), files[2])
    s3 = len(holder)
    print "# - {0} matched stars, {1} unmatched in file 1, {2} unmatched in file 2".format(s1,s2,s3)
    
def create_indices(sorted_data, sub_list, c=0):
    """ Creates a set of indices for the sorted_data, at intervals of sub_list
        c is the column from which to index """
    index = []
    for i in range(len(sorted_data[:,c])):
        if (i % sub_list) == 0:
            index.append((sorted_data[i,c], i))
    return index

if __name__ == "__main__":
    d1 = fi.read_data("stars-10-OLD.txt")
    d2 = fi.read_data("stars-10-NEW.txt")
    sort = Compare(d1, d2, c1=0, c2=0, sub=100)
    match_stars(sort, 5.0*arcsec)
    separate_stars(sort, ["matched-10-2.txt", "unmatched-10-old-2.txt", "unmatched-10-new-2.txt"])
    print "# --- Done"


"""
def get_index(min, max, data):
    base = sc.arange(min, (max+1), 1.0)  # "+1" prevents index overflow during matching
    index, k = [], 0
    for i in range(len(data[:,0])):
        if k == (len(base)-1):  break
        if data[i,0] > base[k]:
            index.append(i)
            k = k + 1
    index.append((len(data[:,0])-1))  # fixes index overflow
    return zip(base, index)

def compare(limit=arcsec):
    # Get Data (l b r)
    data1 = fi.read_data("stars-82.txt")
    data2 = fi.read_data("stars-82-coadd2.txt")
    # Sort data by l then b
    #data1.view('i4,i4,i4').sort(order=['f1'], axis=0)  #rearanges 1st column only!
    data1 = data1[data1[:,0].argsort(),:]
    print "# - Single run data sorted"
    data2 = data2[data2[:,0].argsort(),:]
    print "# - Coadd data sorted"
    # Index Arrays
    if sc.ma.min(data1[:,0]) < sc.ma.min(data2[:,0]):  min = sc.ma.min(data1[:,0])
    else:  min = sc.ma.min(data2[:,0])
    if sc.ma.max(data1[:,0]) > sc.ma.max(data2[:,0]):  max = sc.ma.max(data1[:,0])
    else:  max = sc.ma.max(data2[:,0])
    min, max = int(min), int(max) + 1
    print "# - Indexing from {0} to {1}".format(min, max)
    #index1 = get_index(min, max, data1)
    index2 = get_index(min, max, data2)
    # Do Matching
    holder, match = [], 0
    for i in range(len(data1[:,0])):
        l0, b0 = data1[i,0], data1[i,1]
        # get first index that is higher than l0
        k = 0
        while index2[k][0] < l0:  k = k + 1
        # Use the indices of one step before k to k as the search range
        if (k-1) < 0:  low, high = index2[0][1], index2[k][1]
        else:  low, high = index2[(k-1)][1], index2[k][1]
        # Do the match search over the limited range
        for j in range(low, high):
            l1, b1 = data2[j,0], data2[j,1]
            d = ma.sqrt((l1-l0)*(l1-l0) + (b1-b0)*(b1-b0))
            if d < limit:  match=1;  break 
        if match==0:  holder.append([data1[i,0],data1[i,1],data1[i,2]])
        match=0
        if i % 1000 == 0:
            print "Searched through line {0}, {1} with no matches".format(i, len(holder))
    print "# - Found {0} stars in single run without matches in coadd".format(len(holder))
    return sc.array(holder)

if __name__ == "__main__":
    no_matches = compare(10.0*arcsec)
    fi.write_data(no_matches, "no_matches_82.txt", header="l,b,r; orphan stars in single 82")
    print "# --- Done"

"""
