#!/usr/bin/python

import sys, os
sys.path.insert(1, 'build/lib.macosx-10.10-intel-2.7/')
import comet as C    
from math import log10, exp, factorial, log
from time import time
from scipy.stats import binom
from mutation_data import * 
import json

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
    parser.add_argument('-k', '--gene_set_size', type=int, required=True,
                        help='Gene set size.')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Prefix to output directory.')    
    parser.add_argument('-w', '--weight_func', default='exact', help='Which test to perform.',
                        choices=['binom', 'binom-co', 'exact', 'dendrix'])
    

    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)

    return args


def run( args ):
    # read mutation data
    include = white_and_blacklisting(args.patient_whitelist,
            args.patient_blacklist, args.gene_whitelist, args.gene_blacklist)
    gene2include, patient2include = include

    mutation_data = load_mutation_data_w_cutoff(args.mutation_matrix,
            patient2include, gene2include, args.cutoff)
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data

    print '- %s mutation data: %s genes x %s patients' % (args.name, m, n)

    # exhaustive
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
 


if __name__ == "__main__": run( parse_args() )    
