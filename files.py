import csv
import scipy as sc
import math as m
import numpy as np

'''This program reads and writes files of different types, and is easily callable 

March 17, 2010'''

def read_data(filename, delimiter=None, ignore="---", igval=float("Nan"), splitter=None, 
                skip=0, echo=0):
    """Reads files with custom formating: (maybe make these symbols passable to function?)
    # - ignore line
    $ - important outputs
    ! - warnings, print to screen on read
    splitter - data set delimiter; indicates start of new set of data.  If None,
        assumes a single data set.
    skip - integer, skip this number of lines at top of file
    echo - set to "1" to print lines from the file to the screen
    """
    if (filename[-4:]==".csv" and delimiter==None):  delimiter=","
    readfile = open(filename, "r")
    data_series, data, special = [], [], []
    for line in readfile:
        if skip > 0:  skip=skip-1;  continue
        if (line.strip() == ''):  continue
        if (line[0] == "#"):
            if echo ==1:  print line.strip()
            continue
        if (line[0] == "!"): print line.strip();  continue
        if (line[0] == "$"): special.append( (line.strip()).split(delimiter) ); continue
        if splitter != None:
            if (line[0] == splitter): data_series.append(sc.array(data)); data = []; continue
        # If it gets this far, it's a real line of data to be read in
        readlist = line.strip().rstrip(delimiter)  #removes any trailing delimiters
        readlist = readlist.split(delimiter)
        for j in range(len(readlist)):
            if (readlist[j] == ignore): readlist[j]=igval; continue
            if (readlist[j] == ""): continue #Should do more than just 'continue' here?
            readlist[j] = eval(readlist[j].strip() )
        data.append(readlist)
    readfile.close()
    print '# - File {0} successfully read into data'.format(filename)
    #Get the last of the data into data_series if multiple set are present
    if splitter != None:
        if data != []:  data_series.append(data)
        return data_series
    else:
        return sc.array(data)


def quick_read(filename, delimiter=None):
    """Reads files with custom formating: (maybe make these symbols passable to function?)
    # - ignore line
    ! - warnings, print to screen on read
    """
    print "!!! Function 'quick read' is now being tested!!!"
    if (filename[-4:]==".csv" and delimiter==None):  delimiter=","
    readfile = open(filename, "r")
    data = []
    for line in readfile:
        if (line.strip() == ''):  continue
        if (line[0] == "#"):  continue
        if (line[0] == "!"): print line.strip();  continue
        # If it gets this far, it's a real line of data to be read in
        readlist = line.strip().rstrip(delimiter)  #removes any trailing delimiters
        readlist = readlist.split(delimiter)
        for j in range(len(readlist)):
            readlist[j] = eval(readlist[j].strip() )
        data.append(readlist)
    readfile.close()
    print '# - File {0} successfully read into data'.format(filename)
    return sc.array(data)

def super_write(data_out, fileout='output.txt', delimiter='\t', header=''):
    """Write a column-delimited text file from a 2d data array
    Other delimiters may be used by giving a delimiter as input
    inputing a list of data sets will produce a file readable by 
    using 'splitter' in read_data"""
    print "UNFINISHED FUNCTION"
    return -1
    end = '\n'
    length, width = data_out.shape
    writefile = open(fileout, 'w')
    if (header != ''):
        writefile.write(('#' + header + '\n'))
    for i in range(length):
        writestr = ''
        for j in range(width):
            if j == (width - 1):  suffix = '\n'
            else:  suffix = delimiter
            writestr = writestr + str(data_out[i,j]) + suffix
        writefile.write(writestr)
    writefile.close()
    print '#-Data written as', fileout


def read_data_deprecated(filename, delimiter=None, ignore="---"):
    #Ignores null lines and lines starting with '#'.
    #For reading a text data file into an array; delimiter is the column
    #delimiter - blocks of whitespace by default.  Returns an array of same
    #shape as data set. This is slow for large (100k+ lines) files.  Maybe
    #handle large files soon.  ignore defines the 'no data' symbol.  Maybe
    #read in as lists and then transfer to array? <-- use different function for this.
    readfile = open(filename, "r")
    length, width = 0, 0
    for line in readfile:
    #Maybe include a GNUplot-like array delimiter?
        if (line.strip() == ''): continue
        if (line[0] == "#"): continue  
        length = length + 1
        if width == 0:
            if delimiter == None:
                width = len(line.split())
            else:
                width = len(line.split(delimiter))
    readfile.seek(0)
    data_out = (-1000000.0)*sc.ones((length, width), float)    
    l = 0
    for line in readfile:
        if (line.strip() == ''): continue
        if (line[0] == "#"): continue
        if delimiter == None:
            readlist = line.split()
        else:
            readlist = line.split(delimiter)
        for w in range(width):
            if (readlist[w] == ignore): continue
            if (readlist[w] == ""): continue
            data = eval(readlist[w].strip())
            data_out[l, w] = data
        l = l + 1
    readfile.close()
    print '#-File', filename, 'successfully read into data'
    return sc.array(data_out)

def write_data(data_out, fileout='output.txt', delimiter='\t', header=''):
    #Write a column-delimited text file from a 2d data array
    #Other delimiters may be used by giving a delimiter as input
    end = '\n'
    length, width = data_out.shape
    writefile = open(fileout, 'w')
    if (header != ''):
        writefile.write(('#' + header + '\n'))
    for i in range(length):
        writestr = ''
        for j in range(width):
            if j == (width - 1):  suffix = '\n'
            else:  suffix = delimiter
            writestr = writestr + str(data_out[i,j]) + suffix
        writefile.write(writestr)
    writefile.close()
    print '#-Data written as', fileout

def append_data(data_out, fileout='output.txt', delimiter='\t', note=None):
    #Append a column-delimited text file from a 2d data array to a pre-existing file
    #Other delimiters may be used by giving a delimiter as input
    end = '\n'
    length, width = len(data_out[:,0]), len(data_out[0,:])
    writefile = open(fileout, 'a')
    if (note != None):
        writefile.write(('#' + note + '\n'))
    for i in range(length):
        writestr = ''
        for j in range(width):
            if j == (width - 1):  suffix = '\n'
            else:  suffix = delimiter
            writestr = writestr + str(data_out[i,j]) + suffix
        writefile.write(writestr)
    writefile.close()
    return 1

def get_length(file_in):
    """ Returns the length of a file, excluding empty lines or comments """
    readfile = open(file_in, "r")
    length = 0
    for line in readfile:
        if (line.strip() == ''): continue
        if (line[0] == "#"): continue  
        length = length + 1
    return length

def read_csv(filename):
    #Produces empty lines at end of output data!  <--------------------------<--
    #DEPRECATED
    readfile = open(filename, "rb")
    length = len( readfile.readlines() )
    readfile.seek(0)
    csv_data = csv.reader(readfile)
    i = 0
    for row in csv_data:
        if (row[0].strip() == ''): continue
        if (row[0][0] == "#"): continue
        csv_row = row
        if i == 0:
            width = len(csv_row)
            data_out = (-1.0)*sc.ones((length, width), float)
        for j in range(width):
            data_out[i,j] = eval(csv_row[j])
        i = i + 1
    readfile.close()
    print '#-File', filename, 'successfully read into data'
    return data_out

def write_csv(data_out, filename):
    # Write data to a .csv (comma separated value) file
    #DEPRECATED
    fhandle = open(filename, 'w')
    outfile = csv.writer(fhandle, delimiter=',')
    length, width = data_out.shape
    for i in range(length):
        datarow = []
        for j in range(width):
            datarow.append(str(data_out[i,j]))
        outfile.writerow(datarow)
    fhandle.close()
    return 1

"""Builds a list, where each member is the data from a single run.
Each data set, for n streams:
Likelihood[0], q[1], R0[2], epsilon_1[3], mu_1[4], R_1[5], theta_1[6], phi_1[7], sigma_1[8], ...,
epsilon_n[3+(6*n)], mu_n[4+(6*n)], R_n[5+(6*n)], theta_n[6+(6*n)], phi_n[7+(6*n)], sigma_n[8+(6*n)]"""
def read_boinc_results(filename, num_streams):
    readfile = open(filename, "r")
    #output data list, run list, iteration counter, # parameters per run
    data, run, i, parameters = [], [], 1, (3+(6*num_streams))
    for line in readfile:
        if i > parameters:
            data.append(run)
            run = []
            i = 1
        if (line.strip() == ''): continue
        if (line[0] == "#"): continue
        holder = line.split()
        for param in holder:
            run.append(param.strip(' ,'))
            i=i+1
    data.append(run)
    readfile.close()
    print '#---Boinc separation file', filename, 'with', num_streams, 'streams successfully read.'
    return data

def uncomment(filein, fileout='None'):
    #Strips a file of non-numerical lines and saves it as a different file.
    #Proceed non-numerical info lines with '#' character
    #Automatically creates a logical outfile name if none provided
    if fileout == 'None':
        fileout = filein[0:-4] + '_data' + filein[-4:]
    readfile = open(filein, "r")
    writefile = open(fileout, 'w')
    #length = len( readfile.readlines() )
    #readfile.seek(0)
    for line in readfile:
        #linestr = line #readfile.readline()
        if (line[0] != '#'):
            writefile.write(line)
    readfile.close()
    writefile.close()
    print '-Numerical Data written as', fileout
    return 1


'''I think this function is crap; remove it?'''
def read_line_titles(filename, delimiter='None'):
    #Reads in a data set, where the first column is a string-type title for
    #the data row.
    readfile = open(filename, "r")
    length = len( readfile.readlines() )
    readfile.seek(0)
    readstr = readfile.readline()
    if delimiter == 'None':
        readlist = readstr.split()
    else:
        readlist = readstr.split(delimiter)
    width = len(readlist) 
    readfile.seek(0)
     #'-1' on next line to account for the removal of the line title
    data_out = (-1.0)*sc.ones((length, (width-1)), float)
    title_out = []
    l = 0
    for line in readfile:
        if delimiter == 'None':
            readlist = line.split()
        else:
            readlist = line.split(delimiter)
        for w in range(width):
            if (w == 0):
                title_out.append(readlist[w])
            else:
                data = eval(readlist[w])
                data_out[l, (w-1)] = data
        l = l + 1
    readfile.close()
    print '-File', filename, 'successfully read data and titles, with', \
    len(data_out), 'elements.'
    return title_out, data_out

