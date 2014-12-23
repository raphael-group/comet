
import sys, os
import itertools

# Import c-library
sys.path.insert(1,'build/lib.macosx-10.10-intel-2.7/')
import comet as C
from math import log10, exp, factorial, log
import math
from time import time
import scipy as sp
from mutation_data import * 
import json, re
import convergence as conv
import gc

# Mappings from program arguments to chars passed into C program as flags
EXACT, DENDRIX, BINOM, BINOM_CO = 'e', 'd', 'b', 'c'
weight_function_chars = {
    'exact' : EXACT,
    'dendrix' : DENDRIX,
    'binom' : BINOM,
    'binom-co': BINOM_CO
    }


def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Runs CoMEt to find the optimal set M  '\
                  'of k genes for the weight function \Phi(M).'

    parser = argparse.ArgumentParser(description=description)
    
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-n', '--name', required=True,
                        help='Name for output.')
    # mutation data
    parent_parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parent_parser.add_argument('-c', '--cutoff', type=int, default=0, 
                        help='Minimum gene mutation frequency.')
    parent_parser.add_argument('-p', '--patient_whitelist', default=None,
                        help='File of patients to be included.')
    parent_parser.add_argument('-bp', '--patient_blacklist', default=None,
                        help='File of patients to be excluded.')
    parent_parser.add_argument('-g', '--gene_whitelist', default=None,
                        help='File of genes to be included.')
    parent_parser.add_argument('-bg', '--gene_blacklist', default=None,
                        help='File of genes to be excluded.')
    # output
    parent_parser.add_argument('-o', '--output_dir', required=True,
                        help='Prefix to output directory.')    
    
    subparsers = parser.add_subparsers(help='sub-command help', dest='approach')    
    
    parser_pipeline = subparsers.add_parser('pipeline', parents=[parent_parser], help='pipeline help')
    parser_exhaustive = subparsers.add_parser('exhaustive', parents=[parent_parser], help='exhaustive help')
    
    # Comet parameters
    parser_pipeline.add_argument('-ks', '--gene_set_sizes', nargs="*", type=int, required=True,
                        help='Gene set sizes (length must be t). This or -k must be set. ')
    parser_pipeline.add_argument('-N', '--num_iterations', type=int, default=pow(10, 3),
                        help='Number of iterations of MCMC.')
    parser_pipeline.add_argument('-NStop', '--n_stop', type=int, default=pow(10, 8),
                        help='Number of iterations of MCMC to stop the pipeline.')                  
    parser_pipeline.add_argument('-s', '--step_length', type=int, default=100,
                        help='Number of iterations between samples.')
    parser_pipeline.add_argument('-init', '--initial_soln', nargs="*", 
                        help='Initial solution to use.')
    # new parameters 
    parser_pipeline.add_argument('-acc', '--accelerator', default=1, type=int,
                        help='accelerating factor for target weight')
    parser_pipeline.add_argument('-sub', '--subtype', default=None, help='File with a list of subtype for performing subtype-comet.')
    parser_pipeline.add_argument('-r', '--num_initial', default=1, type=int,
                        help='number of different initial start of MCMC.')
    parser_pipeline.add_argument('--exact_cut', default=0.001, type=float, 
                        help='Maximum accumulated table prob. to stop exact test.')    
    parser_pipeline.add_argument('--binom_cut', type=float, default=0.005,
                        help='Minumum pval cutoff for CoMEt to perform binom test.')
    parser_pipeline.add_argument('-nt', '--nt', default=10, type=int, 
                        help='Maximum co-occurrence cufoff to perform exact test.')
    parser_pipeline.add_argument('-tv', '--total_distance_cutoff', type=float, default=0.005,
                        help='stop condition of convergence (total distance).')
    
    
    # Exhaustive parameters
    parser_exhaustive.add_argument('-k', '--gene_set_size', type=int, required=True,
                        help='Gene set size.')
    parser_exhaustive.add_argument('-w', '--weight_func', default='exact', help='Which test to perform.',
                        choices=['binom', 'binom-co', 'exact', 'dendrix'])
    
    
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

def iter_num (prefix, num_iterations, ks, name):

	if num_iterations >= 1000000000: iterations = "%sB" % (num_iterations / 1000000000)
	elif num_iterations >= 1000000: iterations = "%sM" % (num_iterations / 1000000)
	elif num_iterations >= 1000: iterations = "%sK" % (num_iterations / 1000)
	else: iterations = "%s" % num_iterations
	
	prefix += ".k%s.%s.%s" % ("".join(map(str, ks)), 'pipeline', iterations)
	return prefix

def call_multidendrix(mutation_data, k, t):
    import multi_dendrix as multi_dendrix
    alpha, delta, lmbda = 1.0, 0, 1 # default of multidendrix
    gene_sets_w_weights = multi_dendrix.ILP( mutation_data, t, k, k, alpha, delta, lmbda)	
    multiset = set()
    for gene_set, W in gene_sets_w_weights:
        multiset.update(gene_set)

    return multiset


def printParameters(args, ks, finaltv):
    opts = vars(args)
    opts['total distance'] = finaltv
    prefix = iter_num(args.output_dir + args.name + '.para', args.num_iterations, ks, args.name)
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

def pipeline_run( mutation_data, args):

    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    # read subtype information if args.subtype is specified
    if args.subtype:
        sub_set = [l.rstrip() for l in open(args.subtype)]
        print '- %s subtype information: %s subtypes' % (args.name, len(sub_set))
    else:
        sub_set = list()

    r_c   = args.num_initial 
    t_c   = 1

    t     = len(args.gene_set_sizes) # number of pathways
    ks    = args.gene_set_sizes      # size of each pathway
    N     = args.num_iterations      # number of iteration 
    N_inc = 1.5                 # increamental for non-converged chain
    s     = args.step_length         # step
    N_stop = args.n_stop
    acc = args.accelerator
    nt = args.nt
    hybrid_cutoff = args.binom_cut
    
    # Precompute factorials
    C.precompute_factorials(max(m, n))
    
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
    
    # output CoMEt MCMC results into .tsv
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
                
        open("%s.tsv" % iter_num(args.output_dir + args.name + '.sum', N*(run_num), ks, args.name), "w").write( "\n".join(tbl) )

def exhaustive_run(mutation_data, args):
    # exhaustive
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    iP2G, iG2P, iG2numMuts, g2index, index2g = convert_mutation_data_for_c(genes, patients, mutation2patients,
                                                         patient2mutations)
    k = args.gene_set_size
    pvalthresh = 1.1
    w = args.weight_func
    C.precompute_factorials(max(m, n))
    C.set_weight(weight_function_chars[w])

    res = C.exhaustive(k, m, n, iP2G, iG2numMuts, pvalthresh) 

    solns, weights, tables, probs = res

    res = zip(solns, weights, tables, probs)
    res.sort(key=lambda arr: arr[1], reverse=True) # sort by weight decreasing
    solns   = [ sorted([genes[g] for g in geneset]) for geneset, w, t, p in res]
    weights = [ w for g, w, t, p in res]
    tables  = [ t for g, w, t, p in res]
    probs   = [ p for g, w, t, p in res]


    output_all((solns, weights, tables, probs), args)


def output_all(res, args):
    '''
    Ouputs results of dendrix/dendrix++ to a tsv file.

    Arguments:
    res - a tuple of 5 lists corresponding to gene sets, weights, tables, p-values
    args - the command line arguments to this program
    and frequencies, the results of either the exhaustive or MCMC algorithms

    Outputs the data to a tsv named <outfile_name>

    '''
    # Parse result into shorter var handles
    solns, weights, tables, probs = res

    # Create an output prefix
    prefix = args.output_dir + "/"
    
    prefix += "%s.k%s.%s.exhaust" % (args.name.lower(), args.gene_set_size, args.weight_func)
    
    # Remove sets with a negative weight
    indices = [i for i in range(len(probs)) if probs[i] != -1]
    solns   = [solns[i] for i in indices]
    weights = [weights[i] for i in indices]
    tables  = [tables[i] for i in indices]
    probs   = [probs[i] for i in indices]
 
    # Output only sets, probs, and freqs as TSV
    open("%s.tsv" % prefix, "w").write( "#Gene set\tP-value\tFreq\tWeight\n" +
        "\n".join(["%s\t%s\t%s" % (", ".join(s), p, w)
                   for s, p, w in zip(solns, probs, weights)])
    )
 

def run( args ):
	# read mutation data
    include = white_and_blacklisting(args.patient_whitelist,
            args.patient_blacklist, args.gene_whitelist, args.gene_blacklist)
    gene2include, patient2include = include

    mutation_data = load_mutation_data_w_cutoff(args.mutation_matrix,
            patient2include, gene2include, args.cutoff)
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    print '- %s mutation data: %s genes x %s patients' % (args.name, m, n)

    if args.approach == 'pipeline':
        pipeline_run( mutation_data, args)
    else:
        exhaustive_run( mutation_data, args)
    
    

       
if __name__ == "__main__": run( parse_args() )    
