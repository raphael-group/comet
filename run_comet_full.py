#!/usr/bin/python

# Load required modules
import sys, os, json, re, time, comet as C, multiprocessing as mp, random
from math import exp
import run_comet_simple as RC

def get_parser():
    # Parse arguments
    import argparse
    description = 'Runs CoMEt on permuted matrices.'
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-o', '--output_directory', required=True,
                        help='Output directory.')
    parser.add_argument('--parallel', default=False, action='store_true',
                        help='Use multiprocessing to run a job on each core.')
    parser.add_argument('-np', '--num_permutations', required=True, type=int,
                        help='Number of permuted matrices to use.')
    # Mutation data
    parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-mf', '--min_freq', type=int, default=0,
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-pf', '--patient_file', default=None,
                        help='File of patients to be included (optional).')
    parser.add_argument('-gf', '--gene_file', default=None,
                        help='File of genes to be included (optional).')
    
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
    parser.add_argument('-r', '--num_initial', default=1, type=int,
                        help='Number of different initial starts to use with MCMC.')
    parser.add_argument('-tv', '--total_distance_cutoff', type=float, default=0.005,
                        help='stop condition of convergence (total distance).')
    
    # Parameters for determining the test to be applied in CoMEt
    parser.add_argument('--exact_cut', default=0.001, type=float,
                        help='Maximum accumulated table prob. to stop exact test.')
    parser.add_argument('--binom_cut', type=float, default=0.005,
                        help='Minumum pval cutoff for CoMEt to perform binom test.')
    parser.add_argument('-nt', '--nt', default=10, type=int,
                        help='Maximum co-occurrence cufoff to perform exact test.')
    
    # Files for subtypes/core-events run
    parser.add_argument('-sub', '--subtype', default=None,
                        help='File with a list of subtype for performing subtype-comet.')
    parser.add_argument('-ce', '--core_events', default=None,
                        help='File with a list of core events for performing subtype-comet.')
    
    # Hidden parameters: users can still use these parameters but they won't show in the options    
    # Parameters for marginal probability graph (optional)
    # File mapping genes/events to new names (optional).
    parser.add_argument('-e', '--event_names', default=None, help=argparse.SUPPRESS)
    # File mapping samples to cancer types.
    parser.add_argument('-st', '--sample_types_file', default=None, help=argparse.SUPPRESS)
    # Minimum edge weight for showing in the graph
    parser.add_argument('-mew', '--minimum_edge_weight', type=float, default=0.001,
                            help=argparse.SUPPRESS)
    # Minimum sampling frequency for a gene set to be included.
    parser.add_argument('-msf', '--minimum_sampling_frequency', type=float, default=50,
                            help=argparse.SUPPRESS)
    # Template file (HTML). Change at your own risk.
    parser.add_argument('-tf', '--template_file', default="comet/src/html/template.html",
                            type=str, help=argparse.SUPPRESS)
    # Maximum standard error cutoff to consider a line
    parser.add_argument('-rmse', '--standard_error_cutoff', default=0.01, type=float,
                            help=argparse.SUPPRESS) 
    
    # Input file with lists of pre-run results.
    parser.add_argument('--precomputed_scores', default=None, help=argparse.SUPPRESS)
    # Accelerating factor for target weight
    parser.add_argument('-acc', '--accelerator', default=1, type=int, help=argparse.SUPPRESS)
    # Flag verbose output
    parser.add_argument('-v', '--verbose', default=True, action="store_true",
                        help=argparse.SUPPRESS)
    # Set the seed of the PRNG.
    parser.add_argument('--seed', default=int(time.time()), type=int,
                        help=argparse.SUPPRESS)
    # Edge swapping parameter.
    parser.add_argument('-q', '--Q', type=int, default=100,
                            help=argparse.SUPPRESS)
    # Keep temp files (CoMEt results and permuted matrices).
    parser.add_argument('--keep_temp_files', required=False, action='store_true', default=False,
                        help=argparse.SUPPRESS)
  
    return parser

def runComet(cometArgs):
    return RC.run( RC.get_parser().parse_args(cometArgs) )

def run( args ):

    # Set up the arguments for a general CoMEt run on real data
    realOutputDir = "{}/comet-results".format(args.output_directory)
    realCometArgs = []
    permuteFlags = ["-np", "--parallel", "--keep_temp_files", "-o"]
    for i, arg in enumerate(sys.argv[1:]):
        if arg not in permuteFlags and sys.argv[i] not in permuteFlags:
            realCometArgs.append( arg )
    
    realCometArgs += [ "-o", realOutputDir, "--noviz"]
	# perform simple run without viz first.
    results = runComet(realCometArgs)
    
    # Load mutation data using Multi-Dendrix and output as a temporary file
    realMutations = C.load_mutation_data(args.mutation_matrix, args.patient_file,
                                     args.gene_file, args.min_freq, args.subtype)
    m, n, genes, patients, geneToCases, patientToGenes, subtypes = realMutations

    if args.verbose:
        print '* Mutation data: %s genes x %s patients' % (m, n)

    # Construct bipartite graph from mutation data
    if args.verbose: print "* Creating bipartite graph..."
    G = C.construct_mutation_graph(geneToCases, patientToGenes)
    if args.verbose:
        print '\t- Graph has', len( G.edges() ), 'edges among', len( G.nodes() ), 'nodes.'

    # reset the arguments for a general CoMEt run on permuted matrices
    cometArgs = []
    permuteFlags = ["-np", "--parallel", "--keep_temp_files", "-m", "-o"]
    for i, arg in enumerate(sys.argv[1:]):
        if arg not in permuteFlags and sys.argv[i] not in permuteFlags:
            cometArgs.append( arg )

    cometArgs.append('--noviz')
    # Create a permuted matrix, and then run it through CoMEt
    import tempfile
    arguments = []
    if args.keep_temp_files:
        directory = args.output_directory
    else:
        directory = tempfile.mkdtemp(dir=".", prefix=".tmp")

    # Generate random seeds for each permutation
    random.seed(args.seed)
    seeds = [ random.randint(0, 2**31-1) for _ in range(args.num_permutations) ]

    for i, seed in enumerate(seeds):
        # Print simple progress bar
        sys.stdout.write("* Running CoMEt on permuted matrices... {}/{}\r".format(i+1, args.num_permutations))
        sys.stdout.flush()

        # Create a permuted dataset and save it a temporary file
        mutations = C.permute_mutation_data(G, genes, patients, seed, args.Q)
        _, _, _, _, geneToCases, patientToGenes = mutations
        adj_list = [ p + "\t" + "\t".join( sorted(patientToGenes[p]) ) for p in patients ]

        permutation_file = "{}/permuted-matrix-{}.m2".format(directory, i+1)
        with open(permutation_file, 'w') as outfile: outfile.write('\n'.join(adj_list))

        # Add the new arguments
        permuteArgs = map(str, cometArgs)
        permuteArgs += [ "-m", permutation_file ]
        permuteArgs += [ "-o", "{}/comet-results-on-permutation-{}".format(directory, i+1)]
        arguments.append( permuteArgs )

    if args.parallel:
        pool = mp.Pool(25)
        results = pool.map(runComet, arguments)
        pool.close()
        pool.join()
    else:
        results = [ runComet(permuteArgs) for permuteArgs in arguments ]

    # Find the maximum test statistic on the permuted datasets
    from itertools import islice
    maxStat = 0
    print os.listdir(directory)
    for rf in [ rf for rf in os.listdir(directory) if rf.startswith("comet-results-on-permutation") ]:
        for df in [df for df in os.listdir("{}/{}/results".format(directory, rf)  ) if df.endswith(".tsv")]:            
            with open("{}/{}/results/{}".format(directory, rf, df)) as infile:
                for line in islice(infile, 1, 2):
                    score = float(line.split("\t")[1])        
                    if score > maxStat:
                        maxStat = score

    print "*" * 80
    print "Number of permutations:", args.num_permutations
    print "Max statistic:", maxStat

    # Prepare comet results on real, mutation data, and output directory for viz
    for rf in [rf for rf in os.listdir( "{}/results/".format(realOutputDir) ) if rf.endswith(".tsv")]:
        resultsTable = [l.rstrip() for l in open( "{}/results/{}".format(realOutputDir, rf))]    

    realMutations = (m, n, genes, patients, geneToCases, patientToGenes )
    outputDirViz = realOutputDir + "/viz/"
    C.ensure_dir(outputDirViz)

    # Perform visualization
    C.output_comet_viz(RC.get_parser().parse_args(realCometArgs), realMutations, \
        resultsTable, maxStat, args.num_permutations)

    # Destroy the temporary directory if necessary
    if not args.keep_temp_files:
        import shutil
        shutil.rmtree(directory)


if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )
