#!/usr/bin/python

# Load required modules
import sys, os, json, re, time, comet as C
from math import exp

# Try loading Multi-Dendrix
sys.path.append('third-party/multi-dendrix')
try:
    import multi_dendrix as multi_dendrix    
    importMultidendrix = True
except ImportError:
    importMultidendrix = False
    sys.stderr.write("Note: The Multi-Dendrix Python module could not"\
                     " be found. Using only random initializations...\n")

def get_parser():
    # Parse arguments
    import argparse
    description = 'Runs CoMEt to find the optimal set M  '\
                  'of k genes for the weight function \Phi(M) and '\
                  'visualization webpage. '
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory.')
    # Mutation data
    parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-mf', '--min_freq', type=int, default=0,
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-pf', '--patient_file', default=None,
                        help='File of patients to be included (optional).')
    parser.add_argument('-gf', '--gene_file', default=None,
                        help='File of genes to be included (optional).')
    # CoMEt 
    parser.add_argument('-ks', '--gene_set_sizes', nargs="*", type=int, required=True,
                        help='Gene set sizes (length must be t). This or -k must be set. ')
    # CoMEt MCMC
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


    return parser



def comet(mutations, n, t, ks, numIters, stepLen, initialSoln,
          amp, subt, nt, hybridPvalThreshold, pvalThresh, verbose):
    # Convert mutation data to C-ready format
    if subt: mutations = mutations + (subt, )
    cMutations = C.convert_mutations_to_C_format(*mutations)
    iPatientToGenes, iGeneToCases, geneToNumCases, geneToIndex, indexToGene = cMutations
    initialSolnIndex = [geneToIndex[g] for g in initialSoln]
    solns = C.comet(t, mutations[0], mutations[1], iPatientToGenes, geneToNumCases,
                    ks, numIters, stepLen, amp, nt, hybridPvalThreshold,
                    initialSolnIndex, len(subt), pvalThresh, verbose)

    # Collate the results and sort them descending by sampling frequency
    solnsWithWeights = C.convert_solns( indexToGene, solns )
    def collection_key(collection):
        return " ".join(sorted([",".join(sorted(M)) for M in collection]))

    results = dict()
    # store last soln of sampling for more iterations
    lastSoln = list()
    for gset in solnsWithWeights[-1][0]:
        for g in gset:
            lastSoln.append(g)

    for collection, Ws, Cs in solnsWithWeights:

        key = collection_key(collection)
        if key in results: results[key]["freq"] += 1
        else:
            sets = []
            for i in range(len(collection)):
                M = collection[i]
                W = Ws[i]
                F = Cs[i]
                # extract the probability from the weight,
                # which can also include the accelerator
                P = pow(exp(-W), 1./amp)
                sets.append( dict(genes=M, W=W, num_tbls=F, prob=P) )

            totalWeight  = sum([ S["W"] for S in sets ])
            targetWeight = exp( totalWeight ) if totalWeight < 700 else 1e1000

            results[key] = dict(freq=1, sets=sets, total_weight=totalWeight,
                                target_weight=targetWeight)


    return results, lastSoln


def run( args ):
    ###########################################################################
    # Parse the arguments into shorter variable handles    
    mutationMatrix = args.mutation_matrix
    geneFile = args.gene_file
    patientFile = args.patient_file
    minFreq = args.min_freq
    subtypeFile = args.subtype
    rc    = args.num_initial
    t     = len(args.gene_set_sizes) # number of pathways
    ks    = args.gene_set_sizes      # size of each pathway
    N     = args.num_iterations      # number of iteration
    s     = args.step_length         # step
    NStop = args.n_stop
    acc = args.accelerator
    nt = args.nt
    hybridCutoff = args.binom_cut
    NInc = 1.5                 # increamental for non-converged chain
    tc   = 1

	# Load the mutation data
    mutations = C.load_mutation_data(mutationMatrix, patientFile, geneFile, minFreq, subtypeFile)
    m, n, genes, patients, geneToCases, patientToGenes, subtypes = mutations
    mutations = ( m, n, genes, patients, geneToCases, patientToGenes )


    ###########################################################################
    if args.verbose:
        print 'Mutation data: %s genes x %s patients' % (m, n)

    if args.core_events:
        with open(args.core_events) as f:
            subSet = list( subtypes.union( set( [ l.rstrip() for l in f ] ) ) )
    else:
        subSet = list( subtypes )

    # Precompute factorials
    C.precompute_factorials(max(m, n))
    C.set_random_seed(args.seed)

    # stored the score of pre-computed collections into C
    if args.precomputed_scores:
        C.load_precomputed_scores(args.precomputed_scores, mutations, subSet)

    # num_initial > 1, perform convergence pipeline, otherwise, perform one run only
    if args.num_initial > 1:
        # collect initial soln from users, multidendrix and random.
        initialSolns, totalOut = C.initial_solns_generator(args.num_initial, \
            mutations, ks, args.initial_soln, subSet, \
            importMultidendrix, multi_dendrix)
        runN = N
        while True:
            lastSolns = list()
            for i in range(len(initialSolns)):
                init = initialSolns[i]
                outresults, lastSoln = comet(mutations, n, t, ks, runN, s, \
                    init, acc, subSet, nt, hybridCutoff, args.exact_cut, args.verbose)                
                C.merge_runs(totalOut[i], outresults)
                lastSolns.append(lastSoln)

            finalTv = C.discrete_convergence(totalOut, int(N/s))
            print finalTv, N

            newN = int(N*NInc)
            if newN > NStop or finalTv < args.total_distance_cutoff:
                break
            runN = newN - N
            N = newN
            initialSolns = lastSolns

        runNum = len(totalOut)
        results = C.merge_results(totalOut)
        
    else:
        init = list()
        outresults, lastSoln = comet(mutations, n, t, ks, N, s, \
            init, acc, subSet, nt, hybridCutoff, args.exact_cut, args.verbose)
        results = outresults
        runNum = 1

    C.free_factorials()

    # Output comet results to TSV and website
    collections = sorted(results.keys(), key=lambda S: results[S]["total_weight"], reverse=True)
    C.output_comet(args, mutations, results, collections, ks, N*(runNum), 0, 0)
    
    return [ (S, results[S]["freq"], results[S]["total_weight"]) for S in collections ]

if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )
