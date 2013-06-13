import subprocess as sp
#import glob as glob
import sys

infile = sys.argv[1]
name = infile.split("/")[-1].split(".")[0]

print "\n# --- Beginning compile for file {0} ({1})\n".format(infile, name)
sts1 = sp.call("latex "+infile, shell=True)
print "\n# --- Second pass\n"
sts2 = sp.call("latex "+infile, shell=True)
print "\n# --- Converting to PDF\n"
sts3 = sp.call("dvipdf "+name+".dvi", shell=True)
print "\n# --- Cleaning up extra files\n"
sts4 = sp.call("rm "+name+".log", shell=True)
sts4 = sp.call("rm "+name+".aux", shell=True)
sts4 = sp.call("rm "+name+".dvi", shell=True)
print "\n# --- Compile Complete - output as {0}\n".format(name+".pdf")




"""
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="infile",
                      help="file to be read in, default is first arg", default=sys.argv[1])
    parser.add_option("-t", "--type", dest="plot_type", default='scatter',
                      help="the type of plot to generate.  currently supported: {0}".format(supported))
    parser.add_option("-c", "--columns", dest="columns", default="0:1",
                      help="columns to plot. Inputs depend on plot type: {0}".format(supported_inputs))
    parser.add_option("-x", "--xlabel", dest="xlabel", default="x",
                      help="label for the x axis, accepts latex formatting?")
    parser.add_option("-y", "--ylabel", dest="ylabel", default="y",
                      help="label for the y axis, accepts latex formatting?")
    parser.add_option("--xlim", dest="xlim", default=None,
                      help="range for the x axis - low:high")
    parser.add_option("--ylim", dest="ylim", default=None,
                      help="range for the y axis - low:high")   
    options, extra_args = parser.parse_args(args)
"""
