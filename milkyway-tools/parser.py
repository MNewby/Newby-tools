import sys

'''python script for quickly parsing MW@home results and returning best params.
Usage:  from the command line, type: python parser.py ./folder/population_file
Returns:  Fitness, then a string of all parameters associated with that fitness
Matthew Newby, November 15, 2010'''

in_file = sys.argv[1]
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
print population[best_i][0], population[best_i][1]
#print population[best_i]