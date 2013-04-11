'''python code for producing tables in the aastex deluxe table format
Matthew Newby, March 8, 2012'''

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
            if type(table.data[i][j])==type(""):  entry = float(table.data[i][j])
            else:  entry = table.data[i][j]
            if roundoff != None:
                dataline = dataline + str(round(entry, roundoff)) + " & "
            else:  dataline = dataline + str(entry) + " & "
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