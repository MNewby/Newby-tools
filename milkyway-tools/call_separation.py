import resultsToParameterFile as PF
import separate_streams as ss
#import sdss_visualizers as sdss
import subprocess as sp
import glob as glob
import sys
import getopt

pfolder = "/home/newbym2/milkyway_parameters_sansSgr/"
sfolder = "/home/newbym2/milkyway_stars_sansSgr/"
ofolder = "/home/newbym2/milkyway_out_sansSgr/"
mwfolder = "/home/newbym2/Desktop/milkywayathome_client/bin/"
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
    new_name = MWname.split("_")
    stripe, num_streams = 0, 0
    for name in new_name:  # Get stripe number
        try:
            stripe = int(name)
            if stripe > 8:  break
        except ValueError:  pass
    if stripe == 0:  print "!!!ERROR - Failure to read stripe number"; sys.exit(2)
    for name in new_name:
        if name.count("s") > 0:  
			try:
				num_streams = int(name[0])
			except ValueError:  pass
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
	if stripe==9:  ss="09"
	else:  ss = str(stripe)
	starfiles = glob.glob(sfolder+"*"+ss+"*")
	if len(starfiles) != 1:  
		print "!!!ERROR - {0} starfiles found matching Stripe {1}".format(len(starfiles), stripe)
		if len(starfiles) > 0:
			for s in starfiles:  print s
			sys.exit(2)
	return starfiles[0]


def do_separation(filename, r):
	names, params = parseInFile(filename)
	jobs, coms = [], []
	for i, name in enumerate(names):
		stripe, num_streams = parseString(name)
		print "# --- Processing Stripe {0}".format(stripe)
		p = makeParamFile(stripe, num_streams, params[i], tag)
		s = getStarFile(stripe)
		#s = sfolder+"/stars-"+stripe+"-"+tag+".txt"  #FIX
		if stripe == 9:  ss="09"
		else:  str(stripe)
		o = ofolder+"out-"+ss+"-"+tag+".txt"
		jobs.append([stripe,p,s,o])
		command = mwfolder+"milkyway_separation -i -a "+p+" -s "+s+" -o "+o
		coms.append(command)
		print command
	print "\n# --- Running Commands\n" 
	if r==0:  print "# --- Verification Mode\n"
	else:  print "# --- For Realsies\n"
	for com in coms:  
		print com
		if r==1:  sts = sp.call(com, shell=True)
	print "# --- Done"
		
		
if __name__ == "__main__":
	realsies=0
	resultsfile = sys.argv[1]
	opts, args = getopt.getopt(sys.argv[2:], 'rt:')
	for opt, arg in opts:
		if opt=="-r":  realsies=1
		elif opt=="-t":  tag=arg
	print "# --- Starting runs, using the following settings:"
	print pfolder 
	print sfolder
	print ofolder
	print tfolder
	print tag
	#do_separation(resultsfile, realsies)
	if realsies==1:  
		ss.sep_lbr(sfolder, ofolder, tag)
		ss.separation("/home/newbym2/Desktop/Newby-tools/milkyway-tools/", tag)
		# VISUALIZER
