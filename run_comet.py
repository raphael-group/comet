#!/usr/bin/python

# Load required modules
import sys, os, json, re, comet as C
from math import exp

set2scores = dict()

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
                        choices=C.weight_function_chars.keys())
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


def convert_solns(index2gene, solns):
    new_solns = []
    for arr in solns:
        arr.sort(key=lambda M: M[-2], reverse=True)
        S = tuple( frozenset([index2gene[g] for g in M[:-2] ]) for M in arr )
        W = [ M[-2] for M in arr ]
        F = [ M[-1] for M in arr ]
        new_solns.append( (S, W, F) )

    return new_solns

def comet(mutations, n, t, ks, num_iters, step_len,
          initial_soln, amp, subt, nt, hybrid_pvalthreshold, pvalthresh, verbose):    
    # Convert mutation data to C-ready format
    if subt: mutations = mutations + (subt)
    cMutations = C.convert_mutations_to_C_format(*mutations)
    iPatientToGenes, iGeneToCases, geneToNumCases, geneToIndex, indexToGene = cMutations
    initial_soln_index = [gene2index[g] for g in initial_soln]
    solns = C.comet(t, mutations[0], mutations[1], iPatientToGenes, geneToNumCases,
                    ks, num_iters, step_len, amp, nt, hybrid_pvalthreshold,
                    initial_soln_index, len(subt), pvalthresh, verbose)

    # Collate the results and sort them descending by sampling frequency    
    solns_w_weights = convert_solns( indexToGene, solns )  
    def collection_key(collection):
        return " ".join(sorted([",".join(sorted(M)) for M in collection]))

    results = dict()
    for collection, Ws, Cs in solns_w_weights:
        
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

            total_weight  = sum([ S["W"] for S in sets ])
            target_weight = exp( total_weight ) if total_weight < 700 else 1e1000

            results[key] = dict(freq=1, sets=sets, total_weight=total_weight,
                                target_weight=target_weight)
        
    
    return results

def iter_num (prefix, num_iterations, ks, weight_func, acc):

	if num_iterations >= 1e9: iterations = "%sB" % (num_iterations / 1e9)
	elif num_iterations >= 1e6: iterations = "%sM" % (num_iterations / 1e6)
	elif num_iterations >= 1e3: iterations = "%sK" % (num_iterations / 1e3)
	else: iterations = "%s" % num_iterations
	
	prefix += ".k%s.%s.%s.%s" % ("".join(map(str, ks)), weight_func, iterations, acc)
	return prefix

def call_multidendrix(mutation_data, k, t):
    import multi_dendrix as multi_dendrix
    alpha, delta, lmbda = 1.0, 0, 1 # default of multidendrix
    gene_sets_w_weights = multi_dendrix.ILP( mutation_data, t, k, k, alpha, delta, lmbda)	
    multiset = set()
    for gene_set, W in gene_sets_w_weights:
        multiset.update(gene_set)

    return multiset

def getRanSets(infile):
    
    base_i = 3
    matchObj = re.match( r'.+\.k(\d+)\..+?', infile)
    for l in open(infile):
        if not l.startswith("#"):
            v = l.rstrip().split("\t")
            j = 0
            for i in range(len(matchObj.group(1))):
                set2scores[v[base_i + j]] = float(v[base_i + j + 1])
                #print v[base_i + j], v[base_i + j + 1]
                j += 3

def printParameters(args, ks, finaltv):
    opts = vars(args)
    opts['total distance'] = finaltv
    prefix = iter_num(args.output_prefix + '.para', args.num_iterations, ks, args.weight_func, args.accelerator)
    with open(prefix + '.json', 'w') as outfile:
        json.dump(opts, outfile)


def merge_results(conv_results):
    total = dict()    
    for results in conv_results:
        for key in results.keys():
            if key in total: 
                total[key]["freq"] += results[key]["freq"]
            else:
                total[key] = results[key]

    return total

def run( args ):
    # Parse the arguments into shorter variable hadnles
    mutation_matrix = args.mutation_matrix
    gene_file = args.gene_file
    patient_file = args.patient_file
    min_freq = args.min_freq
    r_c   = args.num_initial 
    t     = len(args.gene_set_sizes) # number of pathways
    ks    = args.gene_set_sizes      # size of each pathway
    w     = args.weight_func         # weight function 
    N     = args.num_iterations      # number of iteration 
    s     = args.step_length         # step
    N_stop = args.n_stop
    acc = args.accelerator
    nt = args.nt
    hybrid_cutoff = args.binom_cut
    N_inc = 1.5                 # increamental for non-converged chain
    t_c   = 1

	# Load the mutation data
    mutations = C.load_mutation_data(mutation_matrix, patient_file, gene_file, min_freq)
    m, n, genes, patients, mutation2patients, patient2mutations = mutations
    
    if args.subtype:
        with open(open(args.subtype)) as f:
            sub_set = [ l.rstrip() for l in f ]
    else:
        sub_set = list()

    if args.verbose:
        print 'Mutation data: %s genes x %s patients' % (m, n)
    
    # Precompute factorials
    C.precompute_factorials(max(m, n))
    
    # stored the score of pre-computed collections
    if args.ran_genesets: 
        getRanSets(args.ran_genesets)

    # num_initial > 1, perform convergence pipeline, otherwise, perform one run only
    if args.num_initial > 1: 
        # create good collection from multidendrix
        pair_init = list()
        if not args.subtype and ks.count(ks[0])==len(ks):        
            md_init = call_multidendrix(mutation_data, ks, t)
            pair_init.append(list(md_init))	       
            t_c = 1    
            r_c = r_c - 1 
        else:
            t_c = 0
			
        while 1:
            conv_results = list()            
            i = 0

            for i in range(t_c):                
                init = pair_init[i]
                outresults = comet(mutations, n, t, ks, N, s, init, acc, sub_set, nt, hybrid_cutoff, args.exact_cut, True)
                conv_results.append(outresults)

            for j in range(i, r_c): # random initials            
                init = list()            
                outresults = comet(mutations, n, t, ks, N, s, init, acc, sub_set, nt, hybrid_cutoff, args.exact_cut, True)
                conv_results.append(outresults)
            
            final_tv = conv.discrete_convergence(conv_results, int(N/s))
            print final_tv, N
                
            newN = int(N*N_inc)
            if newN > N_stop or final_tv < args.total_distance_cutoff: #iterations < 1B and GR factor is < 0.005
                break        
            N = newN
            del conv_results[:]            
        
        run_num = len(conv_results)   
        results = merge_results(conv_results)
        printParameters(args, ks, final_tv) # store and output parameters into .json
        

    else:
        init = list()            
        outresults = comet(mutations, n, t, ks, N, s, init, acc, sub_set, nt, hybrid_cutoff, args.exact_cut, True)
        results = outresults
        run_num = 1
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

    output_file = "%s.tsv" % iter_num(args.output_prefix + '.sum', N*(run_num), ks, w, args.accelerator)
    with open(output_file, "w") as outfile: outfile.write( "\n".join(tbl) )

       
if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )    
