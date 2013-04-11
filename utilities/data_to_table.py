import sys

''' DEPRECATED - SEE "make_deluxe_table.py" 

    Takes a text file and arguments, and returns text output in the format of a 
    asstex 'deluxe table' object.
    
    Inputs:  Space-delimited data file (with text headers), and the number of each
    column to use. (0 is the left-most column; default uses all columns in file)
    Ex:  .../path/to/file 1 3 4 8
Matthew Newby, RPI, 2010 '''

def make_deluxe_table(argv):
    #Handling arguments
    #Checks that first argument is a file name
    if ( (argv[0] == 'help') or (type(argv[0]) != str) ):
        print
        print "---Usage:"
        print "Inputs:    (path/filename) (list of columns to use)"
        print " The second argument is optional; by default all coulumns will be used.\n"
        sys.exit(2)
    #Looks for second argument
    if (len(argv) > 1):
        columns = []
        #checks that second argument is valid; otherwise ignores it.
        try:
            for arg in argv[1:]:
                columns.append(int(arg))
        except ValueError:
            print "!!! ERROR: second argument should be a list of integers separated by spaces"
            print "!!! that define which columns to use, or left blank to imply the use of all columns."+'\n'
            print "!!! Ignoring second argument, using all columns.\n"
            columns = ['all']
    #If no second argument detected, chooses all columns to plot.
    else:
        columns = ['all']
    
    #Reading file
    readfile = open(argv[0], "r")
    checked_header = 0
    contents = []
    for line in readfile:
        if (line.strip() == ''): continue
        if (line[0] == "#"): continue
        readline = line.strip()
        readline = line.split()
        if checked_header == 0:
            try:
                eval(readline[0])
                print "#-No text header found, ignoring header\n"
                header = None
                checked_header = 1
            except NameError:
                print "#-Text header found, read into table header\n"
                header = readline
                checked_header = 1
                continue
        contents.append(readline)  #maybe add something that keep only sig figs if value is a number?
    #Check taht header and data columns are the same length
    if (header != None):
        if ( len(header) != len(contents[0]) ):
            print "!!! ERROR:  Header length not equal to number of data columns"
            print "!!! Ignoring Header...\n"
            header = None
    #get the number of columns for the final table
    if (columns[0] == 'all'):
        width = len(contents[0])
        columns = range(width)
    else:  width = len(columns)
        
    #Create a AASTex deluxetable
    print "#---Beginning deluxe table formatting:\n"
    justify = ''
    for i in range(width):  justify = justify+'r'
    #r, c, or l, determines column justification; maybe make this a third argv?
    print "\\begin{deluxetable}{"+justify+"}"
    print "\\tabletypesize{\scriptsize}"
    print "%\\rotate"  #uncomment this if you want a landscape-style table
    print "\\tablewidth{0pt}"  #Sets table width to the natural width of the page
    print "\\tablecolumns{"+str(width)+"}"
    print "\\tablecaption{ TEXT \label{ KEY }}"  #set short caption and label
    #Produce Table headings:
    print "\\tablehead{"
    if header == None:
        header = []
        for i in range(width):  header.append( ('HEADING '+str(i+1)) )
    colhead = ''
    for i in range(width):
        colhead = colhead+'\colhead{'+header[i]+'} & '
    print colhead[:-2] + '}'
    print "\startdata"
    for i in range(len(contents)):
        dataline = ""
        for j in range(width):
            dataline = dataline + contents[i][columns[j]] + " & "
        print (dataline[:-2] + "\\\\")
    print "\\enddata"
    print "\\end{deluxetable}"
    return 1

if __name__ == "__main__":
    make_deluxe_table(sys.argv[1:])