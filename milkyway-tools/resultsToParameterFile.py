import sys
import getopt

'''python script for quickly parsing MW@home results and building a paramter file out of
the best params.  See below for batch job usage.
Usage:  
Returns:  
Matthew Newby, November 15, 2010'''


def usage():
    print "Usage:"
    print "Arguments:"
    print "# streams:  1, 2, or 3"
    print "original parameter file:  path/name"
    print "call parser?:  -p, to parse a file, -s to use a string of parameters "
    print "EITHER parser file:  path/name"
    print "OR param string in: q R0 epsilon_1 mu_1 R_1 theta_1 phi_1 sigma_1 ... (can have commas)"
    print "Example: 2 ./p-11-2s.txt -p ./results-11-2s-fix.txt"

'''call parser, parameter string in or file to parse, original parameter file in, number of streams?'''

def results_to_parameter_file(argv):
    #Deal with command line arguments
    try: number_streams = int(argv[0])
    except ValueError:
        print "!!!ERROR: # of streams is not an integer!"
        usage()
        sys.exit(2)
    in_file = argv[1]
    if argv[2] == '-p':
        try: fitness, param_string = parse(argv[3])
        except IOError:
            print "!!!ERROR: Unable to read parse file:", argv[3]
            usage()
            sys.exit(2)
        param_list = param_string.split(',')
    elif argv[2] == '-s':
        param_list = argv[3:]
    #Now make the parameters into numbers
    parameters = []
    for param in param_list:
        try: float(param.strip(' ,') )
        except ValueError:
            print "!!!ERROR: unreadable parameter number", param
            usage()
            sys.exit(2)
        if param.count(".") == 0:  param = param + ".0"
        parameters.append(param.strip(' ,') )
    # Check for the right number of paramters
    param_num = 2 + (6*number_streams)
    if param_num != len(parameters):
        print "!!!ERROR: # of input parameters does not match expected # of parameters!"
        usage()
        sys.exit(2)
    #Now we can start modifying the parameter file:
    out_file = in_file[:-4] + '_ps_sanSgr.txt'
    readfile = open(in_file, 'r')
    writefile = open(out_file, 'w')
    stream = 0
    for line in readfile:
        if line[:21] == 'background_parameters':
            front, back = line.split(':')
            back = back.split(',')
            new_line = front+":"+back[0]+', '+parameters[0]+', '+parameters[1]+','+back[3]
        elif line[:14] == 'stream_weight:':
            p_index = 2 + (stream*6)
            new_line = line.split(':')[0]+': '+parameters[p_index]
        elif line[:17] == 'stream_parameters':
            front, back = line.split(':')
            p_index = range( (3+(6*stream)), (8+(6*stream)) )
            back = ', '.join([ parameters[p_index[0]], parameters[p_index[1]],
                             parameters[p_index[2]], parameters[p_index[3]],
                             parameters[p_index[4]] ])
            new_line = front+': '+back
            stream = stream + 1
        elif line[:8] == "convolve":
            new_line = "convolve: 30"
        elif (line[:16] == "r[min,max,steps]") or (line[:5] == "r_cut"):
            front, back = line.split(':')
            back = back.split(',')
            new_line = front+":"+back[0]+','+back[1]+', '+'160'
        elif (line[:17] == "mu[min,max,steps]") or (line[:6] == "mu_cut"):
            front, back = line.split(':')
            back = back.split(',')
            new_line = front+":"+back[0]+','+back[1]+', '+'320'
        elif (line[:17] == "nu[min,max,steps]") or (line[:6] == "nu_cut"):
            front, back = line.split(':')
            back = back.split(',')
            new_line = front+":"+back[0]+','+back[1]+', '+'32'
        else:  new_line = line
        writefile.write((new_line.strip()+'\n'))
    readfile.close()
    writefile.close()
    print "#---File", in_file, "read and closed"
    print "#---File", out_file, "written and closed"

def parse(argv):
    in_file = argv[0]
    #Can uncomment the next line to hard-code inputs - make sure to comment previous line.
    #in_file = './test/validated_population'
    #Open file
    readfile = open(in_file, "r")
    #Assign fitness values and parameters to their own lists
    fitnesses, parameters = [], []
    for line in readfile:
        if (line.strip() == ''): continue
        chop = line.split('fitness')
        if (len(chop) <= 1): continue
        fitnesses.append(float(chop[1][1:-2]))  #fitnesses are converted to floats
        chop = line.split('[')
        numlist = chop[1].split(']')
        parameters.append( numlist[0] )  #parameters are left in string format
    # put the fitnesses and params together:
    population = []
    for i in range(len(fitnesses)):
        population.append([fitnesses[i], parameters[i]])
    # Search fitnesses for best, instead of sorting - which can cause errors if not done carefully.
    # Ignore fitnesses greater than or equalt to 0.0 - they are unphysical
    best_fit, best_i = -1000.0, 1000
    for i in range(len(population)):
        if ( (population[i][0] > best_fit) and (population[i][0] < 0.0) ):
            best_fit, best_i = population[i][0], i
    readfile.close()
    print "#---File", in_file, "read and closed"
    return population[best_i][0], population[best_i][1]
    #print population[best_i]

def batch_job(filename):
    number_of_streams = 3
    stripes = ['09']+range(10,23)
    # Probably should write script to parse files, to make this part quick
    parameter_strings = load_params(filename)
    for i in range(len(stripes)):
        #template_file = "parameters-"+str(stripes[i])+".txt"
        template_file = "./milkyway_parameters/p-"+str(stripes[i])+"-3s-free.txt"
        results_to_parameter_file([3, template_file, "-s"]+parameter_strings[i])
    return 0

def load_params(filename):
    """ loads a set of best-fit parameters from a file """
    holder = []
    readfile = open(filename, "r")
    for line in readfile:
        if line.strip() == "":  continue
        if line[:2] == "ps":  continue
        if line[:2] == "de":  continue
        temp = line.split(',')[2:]    
        for i in range(len(temp)):
            temp[i] = temp[i].strip("[] \n")
        #print temp, "\n"
        holder.append(temp)
    return holder
    

if __name__ == "__main__":
    #results_to_parameter_file(sys.argv[1:])
    batch_job("../sgrOther/Results_ps.txt")  #ALSO CHANGE OUTFILE NAME (line 56)


"""KILL THIS ONLY HERE FOR REFERENCE"""
#<id>214</id> <result><fitness>-3.085803768449465</fitness> <application>stock_win32_gpu: 0.24 double</application>
#<parameters>[0.5806349359602432, 17.777952813089136, -1.6589215177065593, 199.32391451905664, 38.22096391175482, -1.54, -0.049388464632272124, 5.1, -18.318365950027026, 200.34910512566844, 19.728294613893233, -0.06019590689801757, 6.283185307179586, 0.38968648090427405, -20.0, 228.9466772857599, 9.366517140906211, 1.0728682882382794, 5.089967230081957, 5.9569871275727415]</parameters> <metadata>metadata: i: 214</metadata></result>

'''    try: opts, args = getopt.getopt(argv, "?hp:s:d", ["help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-?", "-h", "--help"):
            usage()
            sys.exit(2)
        elif opt == '-p':
            parse_file = arg
        elif opt == '-s':
            param_string = arg  #might cause problems...  maybe just load a file?
            parse_file = ''
    if parse_file != '':
        fitness, param_string = parse(parse_file)
    return 0
'''
