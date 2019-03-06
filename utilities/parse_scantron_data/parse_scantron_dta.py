from __future__ import print_function
import sys
import glob

"""  USAGE INSTRUCTIONS

Matthew Newby, Feb 2, Mar 2, 2019, Temple University"""

# first arg will be filename.  Other args will be special rules
alph = "EDCBA"  #Does the scantron really encode the numbers backwards??? Must test this with an "E" question...
ID_min = 900000000  #this and the next line define the allowed ranges for Student ID numbers
ID_max = 999999999

""" .dta file structure:
1-9    TUID
10     blank
11-21  Last Name
22     blank
23-28  First Name
29     blank
30     Middle Initial
31     blank?
32-36  Special Codes
37-38  Section
39-41  Score (number matching key)
42+    Problems
#/-    End character """

# Add user input safeguards, glob functionality.

""" is the first question always zero?  Is the '#' and '-' a special value? """

def main(args):
    letters, delimiter, fnames, outnames = parse_args(args)
    fcount = 0
    for fn in fnames:
        infile = open(fn, "r")
        lines  = []
        for x in infile:
            lines.append([x[0:9], x[10:21], x[22:28], x[29], x[30:36],
            x[36:38], x[38:41] ] + list(x[41:].strip()[:-1]) )
        infile.close()
        print("File {0} successfully read".format(fn))
        #Now have data in a list:  TUID, Last Name, First Name, MI, codes, section, score, answers
        #switch to letters if needed
        N_probs = len(lines[0][7:])
        if letters == True:
            for i in range(len(lines)):
                for j in range(7, 7+N_probs):  #correct max index?
                    if lines[i][j].strip() == '':  continue
                    lines[i][j] = alph[int(lines[i][j])]
        head, key = create_top(lines)
        error_check(lines)
        simple_out(head, key, lines, outname=outnames[fcount], delim=delimiter)
        fcount = fcount + 1
    print("{0} Files comepleated. Script finished.".format(fcount))

def parse_args(args):
    letters, delim, delimiter = 0, 0, ","
    fnames, outnames = [], []
    if args != []:
        for arg in args:
            if delim == 1:
                delimiter = arg
                delim = 0
            elif "." in arg:
                if len(arg) < 4:
                    print("Invalid filename: {0}".format(arg))
                elif arg[-4:] == ".dta":  #!!! issue if short flags passed before filenames
                    fnames.append(arg)
                elif arg[-4] == ".":
                    outnames.append(arg)
            elif (arg=="letters") or (arg=="alph") or (arg=="-l"):
                letters = 1
            elif (arg=="delimiter") or (arg=="delim") or (arg=="-d"):
                delim = 1
            else:
                print("Invalid argument: {0}".format(arg))
    if fnames == []:
        fnames = glob.glob("*.dta")
        if fnames == []:
            sys.exit("No .dta file in script folder, exiting without action.\n")
    if len(fnames) > len(outnames):
        if (len(fnames) == 1) and (outnames==[]):
            outnames.append("out.csv")
        else:
            for i in range(len(fnames)-len(outnames)):
                outnames.append("out{0}.csv".format(i+1))
    return (letters, delimiter, fnames, outnames)

def format_line(line, delim=","):
    #line is a single line (list) of text, delim is the delimiter to be used.
    temp = ""
    for i in line:
        temp = temp + i.strip() + delim
    temp = temp[:-1] + "\n"
    return temp

def error_check(lines):
    count = 1
    for line in lines:
        errID, errFN, errLN = 0, 0, 0
        if line[0][:3] == "NNN":  continue
        try:
            ID = int(line[0])
            if (ID < ID_min) or (ID > ID_max):
                errID = 1
        except ValueError:
            errID = 1
        if line[1].strip() == '':  eerLN = 1
        if line[2].strip() == '':  eerFN = 1
        if (errID + errFN + errLN) != 0:
            print("Warning: Possible error on line {0}, {1} {2}, ID={3}".format(
            count, line[2], line[1], line[0]), end="" )
            if errFN == 1:  print("; blank first name", end="")
            if errLN == 1:  print("; blank last name", end="")
            if errID == 1:  print("; invalid ID number", end="")
            print()
        count = count + 1

def simple_out(head, key, lines, outname="out.csv", delim=","):
    fout = open(outname, "w")
    if head != None:
        fout.write(format_line(head, delim) )
    if key != None:
        fout.write(format_line(key, delim) )
    for ll in lines:
        if ll[0][0]=="#":  continue
        if ll[0][:3]=="NNN":  continue
        fout.write(format_line(ll, delim) )
    fout.close()
    print("File {0} created.".format(outname) )

def create_top(lines):
    header = ["TUID", "Last Name", "First Name", "MI", "Codes", "Section", "Score"]
    key_line = ["000000000", "ANSWER", "KEY", "",]
    # build header and key
    no_key = 1
    for line in lines:
        if line[0][0:3] == "NNN":
            ans = line[7:]
            N_probs = len(ans)
            key_line.append(line[4]); key_line.append(line[5]); key_line.append(line[6])
            print(line)
            for Q in line[7:]:
                key_line.append(Q)
            no_key = 0
            break
    if no_key == 1:
        print("--No line begins with 'NNN'. Continuing without key line.")
        key_line = None
    for i in range(N_probs):
        header.append("Q{0}".format((i+1) ) )
    return (header, key_line)

if __name__ == "__main__":
    main(sys.argv[1:])
