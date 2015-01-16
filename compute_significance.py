#!/usr/bin/python

# Load required modules
import sys, os, json, math, argparse, time, networkx as nx, tempfile
import comet as C
from collections import defaultdict
import numpy as np
from itertools import combinations as combos
from math import factorial

###############################################################################
# Compute the max weight from a directory of CoMEt results files
def find_max_weight(permutedDir):
    # Identify the CoMEt files (those ending in .tsv)
    from itertools import islice
    permutedFiles = [ f for f in os.listdir( permutedDir ) if f.lower().endswith(".tsv")]
    maxStat = 0.

    # For each file, take the second line which has the top
    # collection by weight, and extract its weight
    for f in permutedFiles:
        with open("{}/{}".format(permutedDir, f)) as infile:
    	    for line in islice(infile, 1, 2):
	            score = float(line.split("\t")[1])
	            if score > maxStat:
	                maxStat = score

    return maxStat, len(permutedFiles)

###############################################################################
# Delta selection functions
# TODO: HSIN-TA: PLEASE ADD COMMENTS TO NEXT THREE FUNCTIONS BELOW

def delta_plot(obj, outfile, boundEdges, deltaPoint, edgeno):
	# TO-DO: HSIN-TA: WHY DOES THIS TAKE IN `obj` instead of `N`, `deltas`, etc.?
	import matplotlib
	matplotlib.use('Agg') # TO-DO: HSIN-TA: WHY?
	import matplotlib.pyplot as plt, seaborn as sns # TO-DO: HSIN-TA: DO WE NEED SEABORN?

	N = obj["N"]
	deltas = obj["deltas"]
	realEdgeDist = obj["edge_dist"]
	plt.rc('text', usetex=True)
	sns.set(style="darkgrid")
	c1, c2, c3, c4 = sns.color_palette("husl", 4)[:4]
	colors = ['b', 'k', 'c', 'y']
	for i in range(2):
		ax = plt.subplot(1, 2, i+1)
		if i == 0 :
			ax.plot(deltas, realEdgeDist, label="Real", c=c4, marker='.')
			ax.set_xlabel("Marginal Probability $p$")
			ax.set_ylabel("No. Edges with Weight $\ge p$")

		if i == 1 :
			ax.plot(deltas, realEdgeDist, label="Real", c='k', marker='.')
			ax.plot(deltaPoint, edgeno, label="Delta", c='r', marker='.')
			ax.axhline(y=boundEdges, xmin=1e-3, xmax=1, c='y')
			ax.set_xlabel("Log (Marginal Probability $p$)")
			ax.set_ylabel("Log (No. Edges with Weight $\ge p$)")
			ax.set_xscale('log')
			ax.set_yscale('log')

			ax.set_xlim([1e-3, 1])
			ax.set_ylim([1, 1000])

		plt.tight_layout()
	plt.savefig(outfile)

def compute_edge_dist(G, deltas):
	weights = [ d['weight'] for u, v, d in G.edges(data=True) ]
	weights.sort()
	N = len(weights)
	numEdges = []
	i = 0
	for j, d in enumerate(deltas):
		M = N if len(numEdges) == 0 else numEdges[-1]
		while i < N and weights[i] < d:
			M -= 1
			i += 1

		numEdges.append( M )

	return np.array(numEdges)

def choose_delta( deltas, realDist, passPoint, stdCutoff):
	"""Find the delta that slope change > 0 on log-log plot"""
	logY = np.log10(np.array(realDist))
	logX = np.log10(np.array(deltas))

	indices = -1
	maxDiff = -1

	startI = 0.
	# reverse order from large edge weight, set start_i as index passing expected edges
	for i in range(len(logX)-1, 0, -1):
		if realDist[i] > passPoint:
			startI = i
			break

	lastSlope = 0.
	i = startI
	lastSlopeI = 0. # for storing the index of elbow point
	while True:
		# define a line based on small standard mean error
		for j in range(2, len(logX)):
			iS = i-j+1
			iE = i+1

			# linear regression
			regression = np.polyfit(logX[iS:iE] , logY[iS:iE], 1)
			# y=ax+b a:point_on_x, b: point_on_y
			pointOnX = regression[1] / (0 - regression[0])
			pointOnY = regression[1]

			# form line 
			line = regression[1] + regression[0] * logX[iS:iE]
			lineX = (logY[iS:iE]-regression[1])/(regression[0])

			# calculate standard mean error
			errX=math.sqrt(sum((lineX-logX[iS:iE])**2)/(iE-iS))
			errY=math.sqrt(sum((line-logY[iS:iE])**2)/(iE-iS))

			if (errX + errY)/2 > stdCutoff or iS == 0:
				break

		if iS == 0 and lastSlope == 0.: # no elbow point
			print "No elbow point found! Try to lower the standard error cutoff with -rmse < ", stdCutoff
			exit(1)

		# first changing point
		if lastSlope == 0.:
			lastSlope = (logY[i] - logY[iS + 1]) / (logX[i] - logX[iS + 1])
			lastSlopeI = iS + 1
		else:
			# slope change > 0
			if lastSlope - ((logY[i] - logY[iS + 1]) / (logX[i] - logX[iS + 1]) )> 0:
				break
			else:
				if iS == 0: # can't find L elbow point after searching all points.
					print "Can't find L (outer) elbow point! Report the last inner elbow point as delta."
					break
				lastSlope = (logY[i] - logY[iS + 1]) / (logX[i] - logX[iS + 1])
				lastSlopeI = iS + 1

		# reset i
		i = iS + 1

	return deltas[lastSlopeI], realDist[lastSlopeI]

###############################################################################
# Functions for input and output

def load_event_names(eventNamesFile, genes):
	# Load a file that maps genes/events/alterations in the mutation to
	# different names for output (presumably shorter names)
	eventNames = dict( (g, g) for g in genes )
	if eventNamesFile:
		with open(eventNamesFile) as f:
			arrs = [ l.rstrip().split() for l in f if not l.startswith("#") ]
			for oldName, newName in arrs:
				eventNames[oldName] = newName
	return eventNames

def create_table_row(arr, eventNames):
	# Convert the given array into a dictionary that can easily
	# be output in JSON and parsed by Javascript
	row = dict()
	row["freq"] = int(arr[0])
	row["totalWeight"] = format(float(arr[1]), 'g')
	row["genesets"] = []
	allGenes = set()
	for i in range(3, len(arr), 3):
		genes = [ eventNames[g] for g in arr[i].split(", ")]
		row["genesets"].append(dict(genes=genes, pval=format(float(arr[i+1]), 'g')))
		allGenes |= set(genes)
	row["allGenes"] = dict((g, True) for g in allGenes)
	return row

# Hsin-Ta: Please comment nCr and construct_mp_graph. What's pass_point, max_rweight,
# etc?
def nCr(n,r):
	return factorial(n) / factorial(r) / factorial(n-r)

def construct_mp_graph(filename, eventNames, minSamplingFreq, maxPWeight):
	tables = []
	geneSets, freqs = [], []
	N = 0
	count = 0
	maxRWeight = 0
	passPoint = 0
	with open(filename) as f:
		for l in f:
			if l.startswith("#") or float(l.split("\t")[1]) < maxPWeight: continue
			if float(l.split("\t")[1]) > maxRWeight:
				maxRWeight = float(l.split("\t")[1])
				passPoint = sum([nCr(len(G.split(", ")),2) for G in l.rstrip().split("\t")[3::3] ])

			count+=1
			arr = l.rstrip().split("\t")
			freq = float(arr[0])

			geneSets.append([G.split(", ") for G in arr[3::3] ])

			freqs.append(freq)

			if freq >= minSamplingFreq:
				tables.append(create_table_row(arr, eventNames))

	# After loading the sets that will go into the graph,
	# compute edges between each pair of genes.
	edges = defaultdict(float)
	N = sum(freqs)
	for Gs, freq in zip(geneSets, freqs):
		for G in Gs:
			for g1, g2 in combos(G, 2):
				if g1 in eventNames and g2 in eventNames:
					edges[frozenset([g1, g2])] += freq/N

	# Create a graph from the edge list
	mpgEdges = [ (list(S)[0], list(S)[1], dict(weight=w)) for S, w in edges.iteritems()  ]
	MPG = nx.Graph()
	MPG.add_edges_from(mpgEdges)

	return MPG, tables, passPoint


def extract_relevant_edges(edges, minEdgeWeight):
	# Sort edges ascending by weight
	edges = [ (u, v, d) for u, v, d in edges if d['weight'] >= minEdgeWeight ]

	# Make a dictionary mapping 
	weightToEdges = defaultdict(list)
	for u, v, d in edges:
		weightToEdges[d['weight']].append( (u, v) )

	return sorted(weightToEdges.keys(), reverse=True), weightToEdges

def gd3_mutation_data(m, n, genes, patients, geneToCases, patientToGenes,
					  genespace, eventNames, sampleToType=None):
	# Sample information (TODO: More elegant solution later)
	if sampleToType:
		typeToSamples = dict( (t, set()) for t in set(sampleToType.values()))
		for sample, ty in sampleToType.iteritems():
			if sample in patients:
				typeToSamples[ty].add( sample )
		typeToSamples = dict( (t, list(s)) for t, s in typeToSamples.iteritems() )
		sampleToTypes = sampleToType
	else:
		typeToSamples = dict(Cancer=patients)
		sampleToTypes = dict((p, "Cancer") for p in patients)
		sampleTypes = ["Cancer"]

	# Create a mutation matrix in GD3 format
	M = {}
	for g in genespace:
		# Determine the type of variant we're looking at
		if "(A)" in g: sym = "amp"
		elif "(D)" in g: sym = "del"
		else: sym = "snv"

		name = eventNames[g]
		if name in M: muts = M[name]
		else: muts = {}
		for p in geneToCases[g]:
			if p in muts: muts[p].append(sym)
			else: muts[p] = [sym]
		M[name] = muts

	return dict(M=M, sampleToTypes=sampleToTypes, typeToSamples=typeToSamples)

def gd3_graph(MPG, eventNames, minEdgeWeight):
	nodes = MPG.nodes()
	edges = [ dict(source=nodes.index(u), target=nodes.index(v), weight=d['weight'])
			  for u, v, d in MPG.edges(data=True) if d['weight'] >= minEdgeWeight ]
	edges.sort(key=lambda d: d['weight'], reverse=True)
	nodes = [ dict(name=eventNames[n]) for n in nodes ]
	weights = sorted([ d['weight'] for d in edges ])
	return dict(nodes=nodes, edges=edges, weights=weights)

###########################################################################$$$$
# Parse command-line arguments and run
def get_parser():
	# Parse arguments
	desc = "Run CoMEt permutation test to create a marginal probability graph"\
		   " in JSON and web form."
	parser = argparse.ArgumentParser(description=desc)

	# CoMEt output
	parser.add_argument('-i', '--input_file', required=True, type=str,
						help='Input file of CoMEt TSV output.')
	parser.add_argument('-pd', '--permuted_dir', default=None, 
						help='Directory with permuted data in tsv.')

	# General
	parser.add_argument('-o', '--output_directory', required=True, type=str,
						help='Path to output directory.')
	parser.add_argument('-v', '--verbose', default=False, action='store_true',
						help='Flag verbose output.')

	# Mutation data
	parser.add_argument('-m', '--mutation_matrix', required=True,
						help='File name for mutation data.')
	parser.add_argument('-mf', '--min_freq', type=int, default=0,
						help='Minimum gene mutation frequency.')
	parser.add_argument('-pf', '--patient_file', default=None,
						help='File of patients to be included.')
	parser.add_argument('-gf', '--gene_file', default=None,
						help='File of genes to be included.')
	parser.add_argument('-e', '--event_names', default=None,
						help='File mapping genes/events to new names (optional).')
	parser.add_argument('-st', '--sample_types_file', default=None,
						 help='File mapping samples to cancer types.')

	# Marginal probability graph
	parser.add_argument('-mew', '--minimum_edge_weight', type=float, default=0.001,
						help='Minimum edge weight.')
	parser.add_argument('-msf', '--minimum_sampling_frequency', type=float, default=50,
					   help='Minimum sampling frequency for a gene set to be included.')
	parser.add_argument('-tf', '--template_file', default="comet/src/html/template.html",
						type=str, help='Template file (HTML). Change at your own risk.')
	parser.add_argument('-rmse', '--standard_error_cutoff', default=0.01, type=float,
					   help='maximum standard error cutoff to consider a line')

	return parser


def run( args ):
	###########################################################################
	# Parse the arguments into shorter variable handles
	mutationMatrix = args.mutation_matrix
	geneFile       = args.gene_file
	patientFile    = args.patient_file
	eventNamesFile = args.event_names
	minFreq        = args.min_freq
	msf            = args.minimum_sampling_frequency
	inputFile      = args.input_file
	permutedDir    = args.permuted_dir
	sec            = args.standard_error_cutoff
	mew            = args.minimum_edge_weight

	# Load the mutation data
	mutations  = C.load_mutation_data(mutationMatrix, patientFile, geneFile, minFreq)
	m, n, genes, patients, geneToCases, patientToGenes = mutations
	eventNames = load_event_names(eventNamesFile, genes)

	###########################################################################
	# Compute max weight from random data if users provide random data.
	# Otherwise, maxPermutedWeight = 0
	if permutedDir:
		maxPermutedWeight, N = calculate_significance(permutedDir)
	else:
		maxPermutedWeight, N = 0, 0

	###########################################################################
	# Construct marginal probability graph from the input CoMEt results file
	if args.verbose: print "* Constructing marginal probability graph..."

	res = construct_mp_graph( inputFile, eventNames,  msf, maxPermutedWeight )
	MPG, tables, passPoint = res
	edges = MPG.edges(data=True)

	if args.verbose: print "\t- Edges:", len(edges)

	# Choose delta (the minimum edge weight in the marginal probability 
	# graph ) using a heuristic approach that selects delta at first elbow
	# with slope change > 0 using linear regression
	if args.verbose: print "* Choosing delta..."

	deltas = sorted(set( d['weight'] for u, v, d in MPG.edges(data=True)))
	realEdgeDist = compute_edge_dist(MPG, deltas)
	deltaPoint, edgeno = choose_delta(deltas, realEdgeDist, passPoint, sec)

	if args.verbose: print "\t- Delta: ", deltaPoint

	###########################################################################
	# Create the web output

	# Dictionary of web output
	obj = {
			"N": N,
			"mm": ['max-derivative'],
			"deltas": deltas,
			"edge_dist": realEdgeDist,
			"collections": {
			   	"max-derivative": {
			   		"more_extreme": 0,
			   		"pval": [0],
				   	"components": list(),
				   	"delta": deltaPoint
				}
			}
	}

	# TO-DO: HSIN-TA: Please comment, I have no idea what this does!
	collections = obj["collections"]
	deltas = [ dict(delta=collections[m]["delta"], pval=collections[m]["pval"], method=m,
			   cdelta=min(obj["deltas"], key=lambda x:abs(x-collections[m]["delta"])))
			   for m in obj["mm"] ]
	plot=None

	# Write the delta plot to file as an SVG, then load
	# it so we can embed it in the web page
	if args.verbose: print "* Plotting delta curve..."
	tmp = tempfile.mktemp(".svg", dir=".", prefix=".tmp")
	delta_plot(obj, tmp, passPoint, deltaPoint, edgeno)

	# Read in the graph, skipping the first four lines that are
	# extraneous header information that will only confuse a webpage
	with open(tmp) as f: plot = "".join(f.readlines()[4:])
	os.unlink(tmp)
	stats = dict(deltas=deltas, plot=plot, N=N)

	# Combine everything to create the D3 data
	if args.verbose: print "* Creating GD3 data..."
	graphData = gd3_graph(MPG, eventNames, mew)
	genesInResults = MPG.nodes()
	sampleToType = None
	if args.sample_types_file:
		with open(args.sample_types_file) as f:
			sampleToType = dict( l.rstrip().split("\t") for l in f )
	mutations = gd3_mutation_data(*mutations, genespace=genesInResults,
								  eventNames=eventNames, sampleToType=sampleToType)

	# Output the results to an HTML file
	if args.verbose: print "* Outputting..."
	htmlOutput = "{}/index.html".format(args.output_directory)
	with open(args.template_file) as template, open(htmlOutput, "w") as outfile:
		jsonData = json.dumps( dict(graph=graphData, mutations=mutations, tables=tables, stats=stats))
		html = template.read()
		html += "\n<script>\nvar data = {};\ndplusViz(data);\n</script>\n".format(jsonData)
		outfile.write( html )

	# Then copy the required JS and CSS files
	import shutil
	shutil.copyfile("comet/src/js/comet-viz.js", "{}/comet-viz.js".format(args.output_directory))
	shutil.copyfile("comet/src/js/mp-graph.js", "{}/mp-graph.js".format(args.output_directory))
	shutil.copyfile("comet/src/js/bower.json", "{}/bower.json".format(args.output_directory))
	shutil.copyfile("comet/src/css/style.css", "{}/style.css".format(args.output_directory))

if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )