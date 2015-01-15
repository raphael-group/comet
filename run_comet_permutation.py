#!/usr/bin/python

# Load required modules
import sys, os, json, re, time, comet as C, multiprocessing as mp
from math import exp
import run_comet as RC

def get_parser():
    # Parse arguments
    import argparse
    description = 'Runs CoMEt on a directory of permuted matrices.'
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-o', '--output_prefix', required=True,
                        help='Output path prefix (TSV format).')
    parser.add_argument('-v', '--verbose', default=True, action="store_true",
                        help='Flag verbose output.')
    parser.add_argument('--seed', default=int(time.time()), type=int,
                        help='Set the seed of the PRNG.')
    parser.add_argument('--parallel', default=False, action='store_true',
                        help='Use multiprocessing to run a job on each core.')

    # Mutation data
    parser.add_argument('-pmd', '--permuted_matrices_directory', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-mf', '--min_freq', type=int, default=0, 
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-pf', '--patient_file', default=None,
                        help='File of patients to be included (optional).')
    parser.add_argument('-gf', '--gene_file', default=None,
                        help='File of genes to be included (optional).')
    parser.add_argument('-np', '--num_permutations', required=True, type=int,
                        help='Number of permuted matrices to use.')

    # Comet
    parser.add_argument('-ks', '--gene_set_sizes', nargs="*", type=int, required=True,
                        help='Gene set sizes (length must be t). This or -k must be set. ')
    parser.add_argument('-N', '--num_iterations', type=int, default=pow(10, 3),
                        help='Number of iterations of MCMC.')
    parser.add_argument('-NStop', '--n_stop', type=int, default=pow(10, 8),
                        help='Number of iterations of MCMC to stop the pipeline.')                  
    parser.add_argument('-s', '--step_length', type=int, default=100,
                        help='Number of iterations between samples.')
    parser.add_argument('-init', '--initial_soln', nargs="*", 
                        help='Initial solution to use.')
    parser.add_argument('-acc', '--accelerator', default=1, type=int,
                        help='accelerating factor for target weight')
    parser.add_argument('-sub', '--subtype', default=None, help='File with a list of subtype for performing subtype-comet.')
    parser.add_argument('-r', '--num_initial', default=1, type=int,
                        help='Number of different initial starts to use with MCMC.')
    parser.add_argument('--exact_cut', default=0.001, type=float, 
                        help='Maximum accumulated table prob. to stop exact test.')    
    parser.add_argument('--binom_cut', type=float, default=0.005,
                        help='Minumum pval cutoff for CoMEt to perform binom test.')
    parser.add_argument('-nt', '--nt', default=10, type=int, 
                        help='Maximum co-occurrence cufoff to perform exact test.')
    parser.add_argument('-tv', '--total_distance_cutoff', type=float, default=0.005,
                        help='stop condition of convergence (total distance).')
    parser.add_argument('--ran_genesets', default=None, 
                        help='input file with lists of pre-run results.')
    return parser

def runComet(cometArgs):
    return RC.run( RC.get_parser().parse_args(cometArgs) )

def run( args ):
    # Set up the arguments for a general CoMEt run 
    cometArgs = []
    permuteFlags = ["-pmd", "-np", "--parallel"]
    for i, arg in enumerate(sys.argv[1:]):
        if arg not in permuteFlags and sys.argv[i] not in permuteFlags:
            cometArgs.append( arg )
    # Identify the permuted matrix files
    filenames = [ f for f in os.listdir(args.permuted_matrices_directory) ]
    filenames = filenames[:min(args.num_permutations, len(filenames))]
    if len(filenames) < args.num_permutations:
        sys.stderr.write("Fewer permuted matrices ({}) than number of permutations requested ({}).\n".format(len(filenames)), args.num_permutations)
    n = len(filenames)

    # Run each matrix through CoMEt
    arguments = []
    for i, f in enumerate(filenames):
        # Print simple progress bar
        sys.stdout.write("* Running CoMEt on permuted matrices... {}/{}\r".format(i+1, n))
        sys.stdout.flush()

        # Add the new arguments
        permuteArgs = map(str, cometArgs)
        permuteArgs += [ "-m", "{}/{}".format(args.permuted_matrices_directory, f)]
        permuteArgs += [ "-o", "{}/{}".format(args.output_prefix, f)]
        arguments.append( permuteArgs )
            
    if args.parallel:
        pool = mp.Pool(25)
        pool.map(runComet, arguments)
        pool.close()
        pool.join()
    else:
        for permuteArgs in arguments:
            runComet(permuteArgs)

if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )