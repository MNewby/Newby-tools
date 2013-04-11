import math as ma
import numpy as np
import scipy as sc
import files as fi
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import functions as func
import re
import glob

'''python script for building a ganged histgram plot for photo metallicity fits
Matthew Newby, Aug 4,2010'''

class HistFit:
    """ """
    def __init__(self, name="EMPTY", params=[0.0,0.0,0.0,0.0,0.0,0.0],
                 errs=[1.0,1.0,1.0,1.0,1.0,1.0]):
        self.name=name
        self.dist=eval(name[-4:].strip("_ "))
        self.params=params
        self.errs=errs
        self.data=sc.zeros(1)
    def insert_data(self, data):
        self.data = data

def unpack_data(resfile, namestr=None):
    """Get data from a custom file, and unpacks the correct entries into HistFits """
    names, parameters, errors = [], [], []
    readfile = open(resfile,"r")
    count, unpack = 1, 0
    for line in readfile:
        if line.strip() == "":  continue
        if line[0] == "#":  continue
        if namestr == None:  unpack = 1
        else:
            if (re.search(namestr, line) != None):  unpack = 1
        # unpack wanted values
        if (unpack == 1) and (count == 1):
            names.append(line.strip())
            count = 2
            continue
        if (unpack == 1) and (count == 2):
            parameters.append(line.strip())
            count = 3
            continue
        if (unpack == 1) and (count == 3):
            errors.append(line.strip())
            count = 1
            unpack = 0
            continue
        # skip unwanted values
        if (unpack == 0) and (count == 1):
            count = 2
            continue
        if (unpack == 0) and (count == 2):
            count = 3
            continue
        if (unpack == 0) and (count == 3):
            count = 1
            unpack = 0
            continue
    readfile.close()
    # if no matching sets were found, exit
    if len(names) == 0:
        print "!!! No entries matching description {0} were found !!!".format(namestr)
        return 0
    fits = []
    #print names, parameters, errors
    for i in range(len(names)):
        params = (parameters[i].strip("[] ")).split()
        for j in range(len(params)):
            params[j] = eval(params[j])
        errs = (errors[i].strip("[] ")).split()
        for j in range(len(errs)):
            errs[j] = eval(errs[j])
        fits.append(HistFit(names[i], params, errs))
    """ Now get the files loaded in """
    files = glob.glob("./An_Project/*"+namestr+"*.txt")
    for i in range(len(fits)):
        lookfor = str(fits[i].dist)
        # The "+1" will cause it to error out if no matching file is found
        for j in range(len(files)+1):
            if re.search(lookfor, files[j]) != None:  break
        fits[i].insert_data(fi.read_data(files[j]))
    return fits

def make_ganged_plot(fits):
    # organize fits
    for fit in fits:
        if fit.name[-3:] == "7.7":  plot1 = fit
        elif fit.name[-4:] == "10.0":  plot2 = fit
        elif fit.name[-4:] == "15.0":  plot3 = fit
        elif fit.name[-4:] == "20.0":  plot4 = fit
        elif fit.name[-4:] == "25.0":  plot5 = fit
        elif fit.name[-4:] == "30.0":  plot6 = fit
        else:  print "! {0} will not be plotted".format(fit.name)
    plots = [plot1, plot4, plot2, plot5, plot3, plot6]
    # initialize plots
    fig = plt.figure()
    plt.subplots_adjust(hspace=0.001, wspace=0.001)
    subs = []
    # do these one at a time
    for i in range(len(plots)):
        """ Setup axes """
        if i == 0:    sp = plt.subplot(3,2,1)
        elif i == 1:  sp = plt.subplot(3,2,2, sharey=subs[0])
        elif i == 2:  sp = plt.subplot(3,2,3, sharex=subs[0])
        elif i == 3:  sp = plt.subplot(3,2,4, sharex=subs[1], sharey=subs[2])
        elif i == 4:  sp = plt.subplot(3,2,5, sharex=subs[0])
        elif i == 5:  sp = plt.subplot(3,2,6, sharex=subs[1], sharey=subs[4])
        subs.append(sp)
        """ plot data and fits """
        x1, x2 = sc.arange(-4.0, 2.0, 0.01), sc.arange(-4.0, 2.0, 0.01)
        y1 = func.gaussian_function(x1, plots[i].params[:3])
        y2 = func.gaussian_function(x2, plots[i].params[3:])
        newx=sc.arange(-4.0, 2.0, 0.01)
        newy=func.double_gaussian(newx, plots[i].params)
        sp.bar(plots[i].data[:,1], plots[i].data[:,0], 0.1, fill=False)
        sp.plot(x1, y1, 'g-')
        sp.plot(x2, y2, 'g-')
        sp.plot(newx, newy, 'r-')
        plt.xlim(-4.0, 2.0)
        plt.ylim(0.0, 140.0)
        plt.xticks(sc.arange(-4.0, 2.0, 1.0), fontsize=8)
        plt.yticks(sc.arange(0.0, 141.0, 20.0), fontsize=8)
        """ Clean up axes and build labels """
        if i == 0:
            plt.setp(sp.get_xticklabels(), visible=False)
            plt.ylabel("counts", fontsize=8)
        elif i == 1:
            plt.setp(sp.get_xticklabels(), visible=False)
            plt.setp(sp.get_yticklabels(), visible=False)
        elif i == 2:
            plt.setp(sp.get_xticklabels(), visible=False)
            plt.ylabel("counts", fontsize=8)
        elif i == 3:
            plt.setp(sp.get_xticklabels(), visible=False)
            plt.setp(sp.get_yticklabels(), visible=False)
        elif i == 4:
            plt.ylabel("counts", fontsize=8)
            plt.xlabel(r"[Fe/H]", fontsize=8)
        elif i == 5:
            plt.xlabel(r"[Fe/H]", fontsize=8)
            plt.setp(sp.get_yticklabels(), visible=False)
        """ insert text """
        if i == 0:
            plt.text(0.0, 100.0, (plots[i].name[-3:]+" kpc"), fontsize=10)
        else:
            plt.text(0.0, 100.0, (plots[i].name[-4:]+" kpc"), fontsize=10)
    plt.show()
    plt.close('all')
    print '#---Done with data plot'

def print_params(fits):
    for fit in fits:
        if fit.name[-3:] == "7.7":  plot1 = fit
        elif fit.name[-4:] == "10.0":  plot2 = fit
        elif fit.name[-4:] == "15.0":  plot3 = fit
        elif fit.name[-4:] == "20.0":  plot4 = fit
        elif fit.name[-4:] == "25.0":  plot5 = fit
        elif fit.name[-4:] == "30.0":  plot6 = fit
        else:  print "! {0} will not be plotted".format(fit.name)
    plots = [plot1, plot2, plot3, plot4, plot5, plot6]
    for plot in plots:
        print plot.name
        print "{0}, {1}, {2}, {3}, {4}, {5}".format(plot.params[0], plot.params[1],
                plot.params[2], plot.params[3], plot.params[4], plot.params[5])
        print "{0}, {1}, {2}, {3}, {4}, {5}".format(plot.errs[0], plot.errs[1],
                plot.errs[2], plot.errs[3], plot.errs[4], plot.errs[5])

if __name__ == "__main__":
    fits = unpack_data("histfit_results.txt", "noF1")
    make_ganged_plot(fits)
    print "# --- Done"