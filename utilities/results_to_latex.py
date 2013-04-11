import math as ma
import numpy as np
import scipy as sc
import matplotlib
import matplotlib.pyplot as plt
#import functions as func
#import astro_coordinates as coor

'''python code for producing tables for the sgrNorth paper.
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

def deluxe_table(table, tofile=None):
    #table = DataTable
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
            dataline = dataline + str(round(table.data[i][j],3)) + " & "
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
        
if __name__ == "__main__":
    import files as fi
    data1 = fi.read_data("../sgrnorth_paper/results_BG_easyread.txt", ",")[:,1:3]
    data2 = fi.read_data("../sgrnorth_paper/results_errors.txt")[:,:2]
    data3 = sc.zeros((len(data1[:,0]),4) )
    data3[:,0], data3[:,1], data3[:,2], data3[:,3] = data1[:,0], data2[:,0], data1[:,1], data2[:,1]
    print data3
    table = DataTable(data3, justify="cccc")
    deluxe_table(table)