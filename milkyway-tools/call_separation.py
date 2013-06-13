import resultsToParameterFile as PF
import subprocess as sp
import glob as glob
import sys
import getopt

pfolder = "/home/newbym2/milkyway_parameters_sansSgr/"
sfolder = "/home/newbym2/milkyway_stars_sansSgr/"
ofolder = "/home/newbym2/milkyway_out_sansSgr/"
mwfolder = "/home/newbym2/milkywayathome_client/bin/"
tfolder = "/home/newbym2/milkyway_parameters_sansSgr/" #template folder
tag=""

def Usage():
	print "Stuff"

def parseInFile(filename):
	"""Reads in a file, formated with the search name in one line,
		then comma-separated parameter strings in a line below that.
		names and parameters should alternate.  If any of this is 
		not true, this function will return an error."""
	names, params = [], []
	results = open(filename, "r")
	for line in results:
		if line.strip()=="":  continue
		if line.strip()[0]=="#":  continue
		if line.strip()[:2]=="de" or line.strip()[:2]=="ps":
			names.append(line.strip().split(",")[0])
			continue
		if (line.count(",") > 0) and (len(names) > len(params)):
			params.append(line.strip())
			continue
		else:  print "!!!ERROR - Invalid results file"; sys.exit(2) 
	results.close()
	if len(names) != len(params):  print "!!!ERROR - Inputs not equal"; sys.exit(2)
	print "Running {0} sets of results:".format(len(params))
	for name in names:  print name
	print
	return names, params

def parseString(MWname):
    """ Parse an MW@H search name string and return the result"""
    new_name = MWname.split("_-")
    for name in new_name:  # Get stripe number
        try:
            stripe = int(name)
            if stripe > 8:  break
        except ValueError:  pass
    for name in new_name:
        if name.count("s") > 0:  num_streams = int(name.strip()[0]);  break
    return stripe, num_streams

def makeParamFile(stripe, num_streams, paramstr, tag=""):
	""" looks for the most valid paramterfile in the parameters directory, 
		then uses it as a template to build a param file for the results."""
	paramfiles = glob.glob(pfolder+"*"+str(stripe)+"*")
	p = None
	if len(paramfiles) > 1:  
		for pfile in paramfiles:  
			if str(num_streams)+"s" in pfile:  
				p = pfile
				if tag in pfile:  break
		if p==None:  p=paramfiles[0]
	elif len(paramfiles)==0:
		print "!!!ERROR - No matching param file for stripe {0}".format(stripe)
		sys.exit(2)
	else: p=paramfiles[0]
	print "Generating parameterfile for stripe {0} using {1}".format(stripe, p.split("\\")[-1])
	PF.results_to_parameter_file([num_streams, p, "-s"]+paramstr.split(",") )
	pname = p[:-4] + '_new.txt'
	return pname

def getStarFile(stripe, tag=""):
	"""Right now, this only works if only one starfile in the target 
		directory matches the stripe number"""
	starfiles = glob.glob(sfolder+"*"+str(stripe)+"*")
	if len(starfiles) != 1:  
		print "!!!ERROR - {0} starfiles found matching Stripe {1}".format(len(starfiles))
		if len(starfiles) > 0:
			for s in starfile:  print s
		sys.exit(2)
	return starfile[0]


def do_separation(filename):
	print "# --- Starting runs, using the following settings:"
	print pfolder 
	print sfolder
	print ofolder
	print tfolder
	print tag
	names, params = parseInFile(filename)
	for i, name in enumerate(names):
		stripe, num_streams = parseString(name)
		print "# --- Starting Stripe {0}".format(stripe)
		p = makeParamFile(stripe, num_streams, params[i], tag)
		s = getStarFile(stripe)
		#s = sfolder+"/stars-"+stripe+"-"+tag+".txt"  #FIX
		o = ofolder+"/out-"+tag+"-"+wedge+".txt"
		command = mwfolder+"/milkyway_separation -i -a "+p+" -s "+s+" -o "+o
		print command
		#sts = sp.call(command, shell=True)
	print "# --- Done"
		
		
if __name__ == "__main__":
	args = sys.argv[1:]
	do_separation(args[0])
