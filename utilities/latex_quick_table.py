import sys

class DataTable:
    """ Contains the headers, data and footnotes for a latex table"""
    def __init__(self, data, justify="r", title="TITLE", label="LABEL", heading=None):
        self.justify = justify
        self.title = title
        self.label = label
        self.heading = heading
        self.data = data
        self.columns = len(data[0])
        self.rows = len(data)

def deluxe_table(table=DataTable, roundoff=None, tofile=None):
    lines = []
    lines.append("\\begin{deluxetable}{"+table.justify+"}")
    lines.append("\\tabletypesize{\scriptsize}")
    lines.append("%\\rotate") #uncomment this if you want a landscape-style table
    lines.append("\\tablewidth{0pt}") #Sets table width to the natural width of the page
    lines.append("\\tablecolumns{"+str(table.columns)+"}")
    lines.append("\\tablecaption{ "+table.title+" \label{ "+table.label+" }}")  #set short caption and label
    #Produce Table headings:
    lines.append("\\tablehead{")
    if table.heading == None:  colhead = "HEADERS }"
    else:
        colhead = ''
        for i in range(table.columns):
            colhead = colhead+'\colhead{'+table.heading[i]+'} & '
        colhead = colhead[:-2] + '}'
    lines.append(colhead)
    lines.append("\\startdata")
    for i in range(table.rows):
        dataline = ""
        for j in range(table.columns):
			entry = table.data[i][j]
			dataline = dataline + round_str(entry, roundoff) + " & "
        dataline = dataline[:-2] + "\\\\"
        lines.append(dataline)
    lines.append("\\enddata")
    lines.append("\\end{deluxetable}")
    # now, output it
    if tofile==None:
        print "% --- START OF LATEX-FORMATTED DELUXE TABLE --- "
        for line in lines:
            print line
        print "% --- END OF LATEX-FORMATTED DELUXE TABLE --- "
    else:
        writefile = open(fileout, 'w')
        for line in lines:
            writefile.write(line)
        writefile.write("\n")
        writefile.close()
        print "# --- Table output to {0}".format(tofile)

def round_str(s, roundoff):
	"""Rounds a number in string format, avoiding pure strings and ints - hopefully"""
	if roundoff==None:  p=s
	elif s.count(".") != 1:  p=s
	else:  p = str(round(float(s), roundoff) )  
	return p

def quick_table(args):
	filein = args[0]
	if len(args) < 2:  delimiter=","
	else:  delimiter = args[1]
	if len(args) < 3:  roundoff = None
	else: roundoff = int(args[2])
	stuff = open(filein, "r")
	holder = []
	while 1:
		line = stuff.readline()
		if line=="\n" or line=="":  
			if holder==[]:  continue
			just = ""
			for i in range(len(head)):  just = just+"c"
			table = DataTable(holder, justify=just, title=title, 
				label="LABEL", heading=head)
			deluxe_table(table, roundoff=roundoff)
			holder = []
			print "\n"
			if line=="":  break
		elif line[0]=="@":  
			hold = line.strip().split(delimiter)
			title = hold[0][1:].strip()
			head = hold[1:]
		elif line[0]=="#":  continue
		else:  holder.append(line.strip().split(delimiter))
	stuff.close()


if __name__ == "__main__":
    args=sys.argv[1:]
    quick_table(args)
    #import files as fi
    #data1 = fi.read_data("../sgrnorth_paper/results_BG_easyread.txt", ",")[:,1:3]
    #data2 = fi.read_data("../sgrnorth_paper/results_errors.txt")[:,:2]
    #data3 = sc.zeros((len(data1[:,0]),4) )
    #data3[:,0], data3[:,1], data3[:,2], data3[:,3] = data1[:,0], data2[:,0], data1[:,1], data2[:,1]
    #print data3
    #table = DataTable(data3, justify="cccc")
    #deluxe_table(table)
