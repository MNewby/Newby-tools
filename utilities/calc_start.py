'''Script to be run with the "PYTHONSTARTUP" environment variable, to turn the
python shell into a calculator without the need to import modules
Matthew Newby, July 11, 2011'''


from math import *
print "# Imported all of math directly to shell"
import numpy as np
print "# Imported numpy as np"
import scipy as sc
print "# imported scipy as sc"
import scipy.constants as con
print "# Imported scipy.constants as con"
#import files as fi
#print "# imported files as fi"
import matplotlib
import matplotlib.pyplot as plt
print "# imported matplotlib.pyplot as plt"
try:
    import astro_coordinates as co
    print "# imported astro_coordinates as co"
except ImportError:
    print "!!! Unable to import astro_coordinates"
print "\n# Entering Interactive Session:"
