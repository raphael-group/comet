#!/usr/bin/python

# Import globally required modules
import sys, os, random, errno
import comet as C

# Parse args
def get_parser():
	# Parse arguments
	import argparse
	description = 'Creates permuted matrices for a given set of CoMEt'\
				  ' mutation data parameters, fixing both the gene '\
				  ' and sample frequencies. Requires NetworkX.'
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-m', '--mutation_matrix', required=True,
						help='File name for mutation data.')
	parser.add_argument('-mf', '--min_freq', type=int, default=0, 
						help='Minimum gene mutation frequency.')
	parser.add_argument('-pf', '--patient_file', default=None,
						help='File of patients to be included.')
	parser.add_argument('-gf', '--gene_file', default=None,
						help='File of genes to be included.')
	parser.add_argument('-o', '--output_dir', required=True,
						help='Name of output directory.')
	parser.add_argument('-s', '--start_index', default=1, type=int,
						help='Start index for name of permuted matrices.')
	parser.add_argument('-n', '--num_matrices', type=int, default=100,
						help='Number of overlaps allowed per pathway.')
	parser.add_argument('-q', '--Q', type=int, default=100,
						help='Edge swapping parameter.')
	parser.add_argument('--verbose', default=False, action='store_true',
						help='Flag verbose mode.')
	return parser

def run(args):
	"""Permutes the given mutation data a given number of times."""
	# Load mutation data using Multi-Dendrix and output as a temporary file
	mutations = C.load_mutation_data(args.mutation_matrix, args.patient_file,
				  				     args.gene_file, args.min_freq)
	m, n, genes, patients, geneToCases, patientToGenes = mutations

	if args.verbose:
		print '* Mutation data: %s genes x %s patients' % (m, n)

	# Make sure output directory exists (credit: http://goo.gl/6J4uX)
	try: os.makedirs(args.output_dir)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(args.output_dir):
			pass
		else: raise
	
	# Construct bipartite graph from mutation data
	if args.verbose: print "* Creating bipartite graph..."
	G = C.construct_mutation_graph(geneToCases, patientToGenes)
	if args.verbose:
		print '\t- Graph has', len( G.edges() ), 'edges among', len( G.nodes() ), 'nodes.'

	# Create permuted matrices and save to file
	matrices = []
	for i in range(args.num_matrices):
		# Simple progress bar
		if args.verbose:
			sys.stdout.write("* Permuting matrices... {}/{}\r".format(i+1, args.num_matrices))
			sys.stdout.flush()
		
		# Permute bipartite graph and output as a patient adjacency list
		mutations = C.permute_mutation_data(G, genes, patients, args.Q)
		_, _, _, _, geneToCases, patientToGenes = mutations
		adj_list = [ p + "\t" + "\t".join( sorted(patientToGenes[p]) ) for p in patients ]
		matrices.append(adj_list)

		filename = "{}/{}.txt".format(args.output_dir, i + args.start_index)
		with open(filename, 'w') as outfile: outfile.write('\n'.join(adj_list))
		
	if args.verbose:
		sys.stdout.write("\r\n")
		sys.stdout.flush()

	return matrices

if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )
