#!/usr/bin/python

# Load required modules
import sys, os, json, re, comet as C
from math import exp

# Try loading Multi-Dendrix
try:
    import multi_dendrix as multi_dendrix
except ImportError:
    importMultidendrix = False
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
    parser.add_argument('-w', '--weight_func', default='exact', help='Which test to perform.',
                        choices=C.weightFunctionChars.keys())
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
    parser.add_argument('--precomputed_scores', default=None, 
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
    initialSolnIndex = [geneToIndex[g] for g in initialSoln]
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

def iter_num (prefix, numIters, ks, weightFunc, acc):

	if numIters >= 1e9: iterations = "%sB" % (numIters / 1e9)
	elif numIters >= 1e6: iterations = "%sM" % (numIters / 1e6)
	elif numIters >= 1e3: iterations = "%sK" % (numIters / 1e3)
	else: iterations = "%s" % numIters
	
	prefix += ".k%s.%s.%s.%s" % ("".join(map(str, ks)), weightFunc, iterations, acc)
	return prefix

def call_multidendrix(mutations, k, t):
    alpha, delta, lmbda = 1.0, 0, 1 # default of multidendrix
    geneSetsWithWeights = multi_dendrix.ILP( mutations, t, k, k, alpha, delta, lmbda)
    multiset = set()
    for geneSet, W in geneSetsWithWeights:
        multiset.update(geneSet)
    return multiset


def initial_solns_generator(r, mutations, ks, assignedInitSoln, subtype):
    runInit = list()
    
    if assignedInitSoln:
        if len(assignedInitSoln) == sum(ks):
            print 'load init soln', "\t".join(assignedInitSoln)
            runInit.append(assignedInitSoln)
        elif len(assignedInitSoln) < sum(ks): # fewer initials than sampling size, randomly pick from genes
            import random
            rand = assignedInitSoln + random.sample(set(mutations[2])-set(assignedInitSoln), sum(ks)-len(assignedInitSoln))
            print 'load init soln with random', rand
        else:
            sys.stderr.write('Too many initial solns for CoMEt.\n')
            exit(1)

    if importMultidendrix and not subtype and ks.count(ks[0])==len(ks): 

        md_init = call_multidendrix(mutations, ks[0], t)
        print ' load multi-dendrix solns', md_init
        runInit.append(list(md_init))

    # assign empty list to runInit as random initials
    for i in range(len(runInit), r):
        runInit.append(list())

    return runInit

def load_precomputed_scores(infile, mutations, subt):
    
    if subt: mutations = mutations + (subt)
    cMutations = C.convert_mutations_to_C_format(*mutations)
    iPatientToGenes, iGeneToCases, geneToNumCases, geneToIndex, indexToGene = cMutations    

    baseI = 3  # sampling freq., total weight, target weight
    setI = 3 # gene set, score, weight function
    
    matchObj = re.match( r'.+\.k(\d+)\..+?', infile)

    loadingT = len(matchObj.group(1)) # determine t:the number of gene sets. 
    for l in open(infile):
        if not l.startswith("#"):
            v = l.rstrip().split("\t")
            j = 0
            for i in range(loadingT):
                gSet = [geneToIndex[g] for g in v[baseI + j].split(", ")]
                C.load_precomputed_scores(float(v[baseI + j + 1]), len(v[baseI + j].split(", ")), int(v[baseI + j + 2]), gSet)
                j += setI

    

def printParameters(args, ks, finaltv):
    opts = vars(args)
    opts['total distance'] = finaltv
    prefix = iter_num(args.output_prefix + '.para', args.num_iterations, ks, args.weight_func, args.accelerator)
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
    w     = args.weight_func         # weight function 
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
    
    # stored the score of pre-computed collections into C
    if args.precomputed_scores: 
        load_precomputed_scores(args.precomputed_scores, mutations, subSet)

    # num_initial > 1, perform convergence pipeline, otherwise, perform one run only
    if args.num_initial > 1:   
        # collect initial soln from users, multidendrix and random.      
        initialSolns = initial_solns_generator(args.num_initial, mutations, ks, args.initial_soln, subSet )
        while True:
            convResults = list()
            for init in initialSolns:
                outresults = comet(mutations, n, t, ks, N, s, init, acc, subSet, nt, hybridCutoff, args.exact_cut, True)  
                convResults.append(outresults)

            finalTv = C.discrete_convergence(convResults, int(N/s))
            print finalTv, N
                
            newN = int(N*NInc)
            if newN > NStop or finalTv < args.total_distance_cutoff: 
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

    outputFile = "%s.tsv" % iter_num(args.output_prefix + '.sum', N*(runNum), ks, w, args.accelerator)
    with open(outputFile, "w") as outfile: outfile.write( "\n".join(tbl) )

if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )    
