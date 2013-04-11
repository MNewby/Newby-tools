import scipy as sc
import sys

'''python script for cutting generic data arrays and returning a subset thereof
Matthew Newby, June 20, 2011'''

def usage():
    print "#--- This script takes in a data array and cut conditions"
    print "data:  any 1 or 2 dimensional numpy (or scipy) array"
    print "cut conditions (args):  series of tuples of the following format:"
    print "\t \t (column#, low_cut, high_cut)"
    print "Where: \n column# is the column to select from"
    print "\t low_cut is the minimum value for this column"
    print "\t high_cut is the maximum value for this column"
    print "Use 'None' for the low/high value to indicate no min/max limit"
    print "\n#-Example: x=cut_data(data, (0, 2.0, 5.0), (2, 1.0, None)) "

def cut_once(data, (col,low_cut,high_cut) ):
    holder=[]
    for i in range(len(data[:,0])):
        if low_cut != None:
            if data[i,col] < low_cut:
                continue
        if high_cut != None:
            if data[i,col] > high_cut:
                continue
        holder.append(data[i,:])
    return sc.array(holder)

""" Main Cutting Function """
def cut_data(data, *args):
    new_data = []
    if len(args) == 0:
        print "!!! NO ARGUMENTS DEFINING CUT !!!\n"
        usage()
        sys.exit(2)
    for arg in args:
        if new_data == []:  new_data = cut_once(data, arg)
        else:  new_data = cut_once(new_data, arg)
    return new_data

if __name__ == "__main__":
    import files as fi
    file_in = "noU_NGC_5466_background.csv"
    data = fi.read_data(file_in, ",")
    data_out = cut_data(data, (2, 20.0, 21.0), (3, None, 20.0))
    fi.write_data(data_out, (file_in[:-4]+"_cut.txt"))
    #print data_out
    #print "length of data:", len(data[:,0]), "; length of cut data:", len(data_out[:,0])
  