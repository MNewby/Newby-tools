import scipy as sc
import math as ma
import numpy as np
import files as fi
import sys
from optparse import OptionParser

"""  This script performs an arbitrary number of cuts on a data set

Matthew Newby, Feb 2, 2012
"""
"""
cutfile = raw_input('Data File:')
column_cut = raw_input('GC distance:')  #column to do cutting in
column = eval(column_cut)
cut_limit_min = raw_input('Min cut:')
cut_min = eval(cut_limit_min)
cut_limit_max = raw_input('Max cut:')
cut_max = eval(cut_limit_max)
"""

def do_one_cut(data, cut_inf):
    col, low, high = cut_inf[0], cut_inf[1], cut_inf[2]
    holder = []
    for i in range(len(data[:,0])):
        if data[i,col] < low: continue
        elif data[i,col] > high:  continue
        else:  holder.append(data[i,:])
    print "# - finished cutting along column {0} from {1} to {2}".format(col, low, high)
    if len(holder)==0:  print "!!! No data exists with those cut values"; sys.exit(2)
    return sc.array(holder)

def do_cuts(infile, delimiter, outfile, cuts):
    # Get data
    data = fi.read_data(infile, delimiter)
    # cut data
    for i in range(len(cuts)):
        data = do_one_cut(data, cuts[i])
    # write data
    fi.write_data(data, outfile, delimiter)
    print "# - Finished with cuts"
    
def parse_options(args):    
    # Parse the data set
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="infile",
                      help="file to be read in, default is first arg", default=sys.argv[1])
    parser.add_option("-d", "--delim", dest="delimiter", default=None,
                      help="the delimiter used in the input file (will also be used in the output file)")
    parser.add_option("-o", "--out", dest="outfile",
                      help="file name for output", default="out.txt")
    parser.add_option("-c", "--cut", action="append", dest='raw_cuts', default=[],
                      help="cuts to make.  format:  [cut_column#]:[low_val]:[high_val]")
    options, extra_args = parser.parse_args(args)
    # Evaluate cut inputs
    cuts = options.raw_cuts
    if len(cuts)==0:  print "!!! invalid number of cuts!"; print cuts; sys.exit(2)
    for i in range(len(cuts)):
        cuts[i] = cuts[i].split(":")
        cuts[i][0] = int(eval(cuts[i][0]))
        cuts[i][1] = eval(cuts[i][1])
        cuts[i][2] = eval(cuts[i][2])
        if cuts[i][1] > cuts[i][2]:
            holder = cuts[i][1]
            cuts[i][1] = cuts[i][2]
            cuts[i][2] = holder
    print cuts
    return options.infile, options.delimiter, options.outfile, cuts

if __name__ == "__main__":
    args = sys.argv[1:]  #parser.parse_args() uses "sys.argv[1:]" by default
    infile, delimiter, outfile, cuts = parse_options(args)
    do_cuts(infile, delimiter, outfile, cuts)
    