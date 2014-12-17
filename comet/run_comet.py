
import sys, os
import itertools

# Import multi-dendrix++
sys.path.insert(1, 'build/lib.macosx-10.10-intel-2.7/')
import comet as C
from math import log10, exp, factorial, log
import math
from time import time
import scipy as sp
from mutation_data import * 
import json, re
import convergence as conv
import gc

set2scores = dict()

def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Runs CoMEt to find the optimal set M  '\
                  'of k genes for the weight function \Phi(M).'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-n', '--name', required=True,
                        help='Name for output.')
    # mutation data
    parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-c', '--cutoff', type=int, default=0, 
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-p', '--patient_whitelist', default=None,
                        help='File of patients to be included.')
    parser.add_argument('-bp', '--patient_blacklist', default=None,
                        help='File of patients to be excluded.')
    parser.add_argument('-g', '--gene_whitelist', default=None,
                        help='File of genes to be included.')
    parser.add_argument('-bg', '--gene_blacklist', default=None,
                        help='File of genes to be excluded.')
    # Comet parameters
    parser.add_argument('-ks', '--gene_set_sizes', nargs="*", type=int, required=True,
                        help='Gene set sizes (length must be t). This or -k must be set. ')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Prefix to output directory.')    
    parser.add_argument('-q', '--quiet', default=True, action="store_true",
                        help='Quiet output flag.')
    parser.add_argument('-N', '--num_iterations', type=int, default=pow(10, 3),
                        help='Number of iterations of MCMC.')
    parser.add_argument('-NStop', '--n_stop', type=int, default=pow(10, 8),
                        help='Number of iterations of MCMC to stop the pipeline.')					
    parser.add_argument('-s', '--step_length', type=int, default=100,
                        help='Number of iterations between samples.')
    parser.add_argument('-w', '--weight_func', default='exact', help='Which test to perform.',
                        choices=['binom', 'exact', 'dendrix'])
    parser.add_argument('-init', '--initial_soln', nargs="*", 
                        help='Initial solution to use.')
    # new parameters 
    parser.add_argument('-acc', '--accelerator', default=1, type=int,
                        help='accelerating factor for target weight')
    parser.add_argument('-sub', '--subtype', default=None, help='File with a list of subtype for performing subtype-comet.')
    parser.add_argument('-r', '--num_initial', default=1, type=int,
                        help='number of different initial start of MCMC.')
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
    
    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)

    return args


def convert_solns(index2gene, solns):
    new_solns = []
    for arr in solns:
        arr.sort(key=lambda M: M[-2], reverse=True)
        S = tuple( frozenset([index2gene[g] for g in M[:-2] ]) for M in arr )
        W = [ M[-2] for M in arr ]
        F = [ M[-1] for M in arr ]
        new_solns.append( (S, W, F) )

    return new_solns

def comet(mutation_data, n, t, ks, num_iters, step_len,
                        initial_soln, amp, subt, nt, hybrid_pvalthreshold, pvalthresh, verbose):    
    # Convert mutation data to C-ready format
    iP2G, iG2P, iG2numMuts, gene2index, index2gene = convert_mutation_data_for_c_with_subtype(mutation_data[2], mutation_data[3], mutation_data[4], mutation_data[5], subt)
    initial_soln_index = [gene2index[g] for g in initial_soln]
    solns = C.comet(t, mutation_data[0], mutation_data[1], iP2G, iG2numMuts, ks, num_iters, step_len, 
                            amp, nt, hybrid_pvalthreshold, initial_soln_index, len(subt), pvalthresh, verbose)

    # Collate the results and sort them descending by sampling frequency    
    solns_w_weights = convert_solns( index2gene, solns )  
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

def iter_num (prefix, num_iterations, ks, weight_func, name, acc):

	if num_iterations >= 1000000000: iterations = "%sB" % (num_iterations / 1000000000)
	elif num_iterations >= 1000000: iterations = "%sM" % (num_iterations / 1000000)
	elif num_iterations >= 1000: iterations = "%sK" % (num_iterations / 1000)
	else: iterations = "%s" % num_iterations
	
	prefix += ".k%s.%s.%s.%s.%s" % ("".join(map(str, ks)), weight_func, iterations, acc)
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
    prefix = iter_num(args.output_dir + args.name + '.para', args.num_iterations, ks, args.weight_func, args.name, args.accelerator)
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
	# read mutation data
    include = white_and_blacklisting(args.patient_whitelist,
            args.patient_blacklist, args.gene_whitelist, args.gene_blacklist)
    gene2include, patient2include = include

    mutation_data = load_mutation_data_w_cutoff(args.mutation_matrix,
            patient2include, gene2include, args.cutoff)
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    
    if args.subtype:
        sub_set = [l.rstrip() for l in open(args.subtype)]
    else:
        sub_set = list()

    
    print '- %s mutation data: %s genes x %s patients' % (args.name, m, n)

    r_c   = args.num_initial 
    t_c   = 1

    t     = len(args.gene_set_sizes) # number of pathways
    ks    = args.gene_set_sizes      # size of each pathway
    w     = args.weight_func         # weight function 
    N     = args.num_iterations      # number of iteration 
    N_inc = 1.5                 # increamental for non-converged chain
    s     = args.step_length         # step
    N_stop = args.n_stop
    acc = args.accelerator
    nt = args.nt
    hybrid_cutoff = args.binom_cut
    
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
                outresults = comet(mutation_data, n, t, ks, N, s, init, acc, sub_set, nt, hybrid_cutoff, args.exact_cut, True)
                conv_results.append(outresults)

            for j in range(i, r_c): # random initials            
                init = list()            
                outresults = comet(mutation_data, n, t, ks, N, s, init, acc, sub_set, nt, hybrid_cutoff, args.exact_cut, True)
                conv_results.append(outresults)
            
            final_tv = conv.discrete_convergence(conv_results, int(N/s))
            print final_tv, N
                
            newN = int(N*N_inc)
            if newN > N_stop or final_tv < args.total_distance_cutoff: #iterations < 1B and GR factor is < 0.005
                break        
            N = newN
            del conv_results[:]            
            #print gc.collect()
            #print gc.garbage

        
        run_num = len(conv_results)   
        results = merge_results(conv_results)
        printParameters(args, ks, final_tv) # store and output parameters into .json
        

    else:
        init = list()            
        outresults = comet(mutation_data, n, t, ks, N, s, init, acc, sub_set, nt, hybrid_cutoff, args.exact_cut, True)
        results = outresults
        run_num = 1
        printParameters(args, ks, 1) 

    C.free_factorials()
    
    # output Comet results into .tsv
    collections = sorted(results.keys(), key=lambda S: results[S]["total_weight"], reverse=True)    		    
    if args.output_dir:   
        # Output only sets, probs, and freqs as TSV
        header = "#Freq\tTotal Weight\tTarget Weight\t"
        header += "\t".join(["Gene set %s (k=%s)\tProb %s\tWeight function %s" % (i, ks[i-1], i, i) for i in range(1, len(ks)+1)])
        tbl = [header]
        for S in collections:
            data = results[S]
            row = [ data["freq"], data["total_weight"], format(data["target_weight"], 'g') ]
            for d in sorted(data["sets"], key=lambda d: d["W"]):
                row += [", ".join(sorted(d["genes"])), d["prob"], d["num_tbls"] ]
            tbl.append("\t".join(map(str, row)))
                
        open("%s.tsv" % iter_num(args.output_dir + args.name + '.sum', N*(run_num), ks, w, args.name, args.accelerator), "w").write( "\n".join(tbl) )

       
if __name__ == "__main__": run( parse_args() )    
