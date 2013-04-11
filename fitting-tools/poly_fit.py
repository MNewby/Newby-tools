import math as m
import numpy as np
import scipy as sc
import files as f
import sys

'''python script for quickly fitting polynomials to data.
See Chapter 7.2, Data Reduction... , Bevington, 2nd Ed.
Matthew Newby, RPI, May 26, 2011'''


def usage():
    print "#- Poly_fit inputs: \n name of data file (1 str), x,y,sigma column #s\
    (3 int), order of polynomial (1 int).\n Can use 'p' for sigma column to get \
    poisson counting errors."
    sys.exit(2)
    return 0

def poly (x_in, N):
    #returns x_in to the Nth power
    holder = 1.0
    for i in range(N):
        holder = holder*x_in
    return holder

def poly_fit(x, y, sig, order, verbose=1):
    n_params = order + 1
    #initialize matrices as arrays
    beta = sc.zeros(n_params, float)
    solution = sc.zeros(n_params, float)
    alpha = sc.zeros( (n_params,n_params), float)
    #Fill Matrices
    for k in range(n_params):
        # Fill beta
        for i in range(len(x)):
            holder = ( y[i]*poly(x[i], k) ) / (sig[i]*sig[i])
            beta[k] = beta[k] + holder
        # Fill alpha
        for l in range(n_params):
            for i in range(len(x)):
                holder = (poly(x[i],l)*poly(x[i],k)) / (sig[i]*sig[i])
                alpha[l,k] = alpha[l,k] + holder
    # Re-type matrices
    beta_m = sc.matrix(beta)
    alpha_m = sc.matrix(alpha)
    #Invert alpha,, then multiply beta on the right by the new matrix
    #epsilon_m = alpha_m.I
    a_m = beta_m * alpha_m.I
    if verbose==1:
        print "beta:\n", beta_m
        print "alpha:\n", alpha_m
        print "best-fit parameter matrix:\n", a_m
    return sc.array(a_m)[0,:]
    
if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) != 5:  usage()
    if (args[0] == "--help") or (args[1] == "help"):  usage()
    import plot_data_function as pdf
    import functions as func    
    data = f.read_data(args[0])
    x, y = data[:,int(args[1])], data[:,int(args[2])]
    if args[3] == 'p':
        sig = func.poisson_errors(y)
        print "#- Using poisson counting errors"
    else:  sig = data[:,int(args[3])]
    results = poly_fit(x,y,sig,int(args[4]))
    pdf.plot_function(x, y, func.Nth_order_polynomial, results, y_err=sig,
                  title=['', 'x', 'y'], save=0, save_name='func_plot')
    print "#-Done"