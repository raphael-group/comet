#!/usr/bin/python

# Load required modules
import sys, os, json, re, time, comet as C, multiprocessing as mp, random
from math import exp
import run_comet as RC

def get_parser():
    # Parse arguments
    import argparse
    description = 'Runs CoMEt on permuted matrices.'
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-o', '--output_directory', required=True,
                        help='Output directory.')
    parser.add_argument('-v', '--verbose', default=True, action="store_true",
                        help='Flag verbose output.')
    parser.add_argument('--seed', default=int(time.time()), type=int,
                        help='Set the seed of the PRNG.')
    parser.add_argument('--parallel', default=False, action='store_true',
                        help='Use multiprocessing to run a job on each core.')
    parser.add_argument('-np', '--num_permutations', required=True, type=int,
                        help='Number of permuted matrices to use.')
    parser.add_argument('--keep_temp_files', required=False, action='store_true', default=False,
                        help='Keep temp files (CoMEt results and permuted matrices).')

    # Mutation data
    parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-mf', '--min_freq', type=int, default=0,
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-pf', '--patient_file', default=None,
                        help='File of patients to be included (optional).')
    parser.add_argument('-gf', '--gene_file', default=None,
                        help='File of genes to be included (optional).')
    parser.add_argument('-q', '--Q', type=int, default=100,
                            help='Edge swapping parameter.')

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
    # Seed Python PRNG
    random.seed(args.seed)

    # Load mutation data using Multi-Dendrix and output as a temporary file
    mutations = C.load_mutation_data(args.mutation_matrix, args.patient_file,
                                     args.gene_file, args.min_freq)
    m, n, genes, patients, geneToCases, patientToGenes = mutations

    if args.verbose:
        print '* Mutation data: %s genes x %s patients' % (m, n)

    # Construct bipartite graph from mutation data
    if args.verbose: print "* Creating bipartite graph..."
    G = C.construct_mutation_graph(geneToCases, patientToGenes)
    if args.verbose:
        print '\t- Graph has', len( G.edges() ), 'edges among', len( G.nodes() ), 'nodes.'

    # Set up the arguments for a general CoMEt run
    cometArgs = []
    permuteFlags = ["-np", "--parallel", "--keep_temp_files", "-m", "-o"]
    for i, arg in enumerate(sys.argv[1:]):
        if arg not in permuteFlags and sys.argv[i] not in permuteFlags:
            cometArgs.append( arg )

    # Create a permuted matrix, and then run it through CoMEt
    import tempfile
    arguments = []
    if args.keep_temp_files:
        directory = args.output_directory
    else:
        directory = tempfile.mkdtemp(dir=".", prefix=".tmp")

    for i in range(args.num_permutations):
        # Print simple progress bar
        sys.stdout.write("* Running CoMEt on permuted matrices... {}/{}\r".format(i+1, args.num_permutations))
        sys.stdout.flush()

        # Create a permuted dataset and save it a temporary file
        seed = random.randint(0, 2**31-1) # generate a new random seed
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
    for rf in [ rf for rf in os.listdir(directory) if rf.startswith("comet-results") ]:
        with open("{}/{}".format(directory, rf)) as infile:
    	    for line in islice(infile, 1, 2):
                score = float(line.split("\t")[1])
                if score > maxStat:
                    maxStat = score

    print "*" * 80
    print "Number of permutations:", args.num_permutations
    print "Max statistic:", maxStat

    # Output the results to files
    with open("{}/comet-stats.json".format(args.output_directory), "w") as outfile:
        output = dict(maxPermutedWeight=maxStat,
                      numPermutations=args.num_permutations,
                      keepTempFiles=args.keep_temp_files,
                      mutationNatrix=args.mutation_matrix,
                      geneFile=args.gene_file, patientFile=args.patient_file,
                      minFreq=args.min_freq, Q=args.Q)
        json.dump( output, outfile, sort_keys=True, indent=4)

    # Destroy the temporary directory if necessary
    if not args.keep_temp_files:
        import shutil
        shutil.rmtree(directory)


if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )
