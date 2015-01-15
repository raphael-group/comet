#!/usr/bin/python

# Load required modules
import sys, os, json, re, time, comet as C
from math import exp

# Try loading Multi-Dendrix
try:
    import multi_dendrix as multi_dendrix
except ImportError:
    sys.stderr.write("Warning: The Multi-Dendrix Python module could not"\
                     " be found. Using only random initializations...\n")

def get_parser():
    # Parse arguments
    import argparse
    description = 'Runs CoMEt to find the optimal set M  '\
                  'of k genes for the weight function \Phi(M).'
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-o', '--output_prefix', required=True,
                        help='Output path prefix (TSV format).')
    parser.add_argument('-v', '--verbose', default=True, action="store_true",
                        help='Flag verbose output.')
    parser.add_argument('--seed', default=int(time.time()), type=int,
                        help='Set the seed of the PRNG.')

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


def convert_solns(indexToGene, solns):
    newSolns = []
    for arr in solns:
        arr.sort(key=lambda M: M[-2], reverse=True)
        S = tuple( frozenset([indexToGene[g] for g in M[:-2] ]) for M in arr )
        W = [ M[-2] for M in arr ]
        F = [ M[-1] for M in arr ]
        newSolns.append( (S, W, F) )

    return newSolns

def comet(mutations, n, t, ks, numIters, stepLen, initialSoln, 
          amp, subt, nt, hybridPvalThreshold, pvalThresh, verbose):    
    # Convert mutation data to C-ready format
    if subt: mutations = mutations + (subt)
    cMutations = C.convert_mutations_to_C_format(*mutations)
    iPatientToGenes, iGeneToCases, geneToNumCases, geneToIndex, indexToGene = cMutations
    initialSolnIndex = [gene2index[g] for g in initialSoln]
    solns = C.comet(t, mutations[0], mutations[1], iPatientToGenes, geneToNumCases,
                    ks, numIters, stepLen, amp, nt, hybridPvalThreshold,
                    initialSolnIndex, len(subt), pvalThresh, verbose)

    # Collate the results and sort them descending by sampling frequency    
    solnsWithWeights = convert_solns( indexToGene, solns )  
    def collection_key(collection):
        return " ".join(sorted([",".join(sorted(M)) for M in collection]))

    results = dict()
    for collection, Ws, Cs in solnsWithWeights:
        
        key = collection_key(collection)
        if key in results: results[key]["freq"] += 1
        else:
            sets = []
            for i in range(len(collection)):
                M = collection[i]
                W = Ws[i]
                F = Cs[i]
                P = exp(-W)
                sets.append( dict(genes=M, W=W, num_tbls=F, prob=P) )

            totalWeight  = sum([ S["W"] for S in sets ])
            targetWeight = exp( totalWeight ) if totalWeight < 700 else 1e1000

            results[key] = dict(freq=1, sets=sets, total_weight=totalWeight,
                                target_weight=targetWeight)
        
    
    return results

def iter_num (prefix, numIters, ks, acc):

	if numIters >= 1e9: iterations = "%sB" % (numIters / 1e9)
	elif numIters >= 1e6: iterations = "%sM" % (numIters / 1e6)
	elif numIters >= 1e3: iterations = "%sK" % (numIters / 1e3)
	else: iterations = "%s" % numIters
	
	prefix += ".k%s.%s.%s" % ("".join(map(str, ks)), iterations, acc)
	return prefix

def call_multidendrix(mutations, k, t):
    alpha, delta, lmbda = 1.0, 0, 1 # default of multidendrix
    geneSetsWithWeights = multi_dendrix.ILP( mutations, t, k, k, alpha, delta, lmbda)
    multiset = set()
    for geneSet, W in geneSetsWithWeights:
        multiset.update(geneSet)
    return multiset

def getRanSets(infile):
    baseI = 3
    matchObj = re.match( r'.+\.k(\d+)\..+?', infile)
    for l in open(infile):
        if not l.startswith("#"):
            v = l.rstrip().split("\t")
            j = 0
            for i in range(len(matchObj.group(1))):
                setToScores[v[baseI + j]] = float(v[baseI + j + 1])
                #print v[base_i + j], v[base_i + j + 1]
                j += 3

def printParameters(args, ks, finaltv):
    opts = vars(args)
    opts['total distance'] = finaltv
    prefix = iter_num(args.output_prefix + '.para', args.num_iterations, ks, args.accelerator)
    with open(prefix + '.json', 'w') as outfile:
        json.dump(opts, outfile)


def merge_results(convResults):
    total = dict()    
    for results in convResults:
        for key in results.keys():
            if key in total: 
                total[key]["freq"] += results[key]["freq"]
            else:
                total[key] = results[key]

    return total

def run( args ):
    # Parse the arguments into shorter variable handles
    mutationMatrix = args.mutation_matrix
    geneFile = args.gene_file
    patientFile = args.patient_file
    minFreq = args.min_freq
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
    mutations = C.load_mutation_data(mutationMatrix, patientFile, geneFile, minFreq)
    m, n, genes, patients, geneToCases, patientToGenes = mutations
    
    if args.subtype:
        with open(open(args.subtype)) as f:
            subSet = [ l.rstrip() for l in f ]
    else:
        subSet = list()

    if args.verbose:
        print 'Mutation data: %s genes x %s patients' % (m, n)
    
    # Precompute factorials
    C.precompute_factorials(max(m, n))
    C.set_random_seed(args.seed)
    
    # stored the score of pre-computed collections
    if args.ran_genesets: 
        getRanSets(args.ran_genesets)

    # num_initial > 1, perform convergence pipeline, otherwise, perform one run only
    if args.num_initial > 1: 
        # create good collection from multidendrix
        pair_init = list()
        if not args.subtype and ks.count(ks[0])==len(ks):        
            mdInit = call_multidendrix(mutations, ks, t)
            pair_init.append(list(mdInit))	       
            tc = 1    
            rc = rc - 1 
        else:
            tc = 0
			
        while 1:
            convResults = list()            
            i = 0

            for i in range(tc):                
                init = pair_init[i]
                outresults = comet(mutations, n, t, ks, N, s, init, acc, sub_set, nt, hybridCutoff, args.exact_cut, True)
                conv_results.append(outresults)

            for j in range(i, rc): # random initials            
                init = list()            
                outresults = comet(mutations, n, t, ks, N, s, init, acc, sub_set, nt, hybridCutoff, args.exact_cut, True)
                conv_results.append(outresults)
            
            final_tv = conv.discrete_convergence(conv_results, int(N/s))
            print final_tv, N
                
            newN = int(N*NInc)
            if newN > NStop or finalTv < args.total_distance_cutoff: #iterations < 1B and GR factor is < 0.005
                break        
            N = newN
            del convResults[:]            
        
        runNum = len(convResults)   
        results = merge_results(convResults)
        printParameters(args, ks, finalTv) # store and output parameters into .json
        

    else:
        init = list()            
        outresults = comet(mutations, n, t, ks, N, s, init, acc, subSet, nt, hybridCutoff, args.exact_cut, True)
        results = outresults
        runNum = 1
        printParameters(args, ks, 1) 

    C.free_factorials()
    
    # Output Comet results to TSV
    collections = sorted(results.keys(), key=lambda S: results[S]["total_weight"], reverse=True)    		    
    header = "#Freq\tTotal Weight\tTarget Weight\t"
    header += "\t".join(["Gene set %s (k=%s)\tProb %s\tWeight function %s" % (i, ks[i-1], i, i) for i in range(1, len(ks)+1)])
    tbl = [header]
    for S in collections:
        data = results[S]
        row = [ data["freq"], data["total_weight"], format(data["target_weight"], 'g') ]
        for d in sorted(data["sets"], key=lambda d: d["W"]):
            row += [", ".join(sorted(d["genes"])), d["prob"], d["num_tbls"] ]
        tbl.append("\t".join(map(str, row)))

    outputFile = "%s.tsv" % iter_num(args.output_prefix + '.sum', N*(runNum), ks, args.accelerator)
    with open(outputFile, "w") as outfile: outfile.write( "\n".join(tbl) )

    return [ (S, results[S]["freq"], results[S]["total_weight"]) for S in collections ]

if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )    
