# output comet results and website
import sys, os, json, math, argparse, time, networkx as nx, tempfile
from collections import defaultdict
import numpy as np
from itertools import combinations as combos
from math import factorial

def ensure_dir(path):
	try: 
		os.makedirs(path)
	except OSError:
		if not os.path.isdir(path):
			raise

def iter_num (prefix, numIters, ks, acc):

	if numIters >= 1e9: iterations = "%sB" % (numIters / int(1e9))
	elif numIters >= 1e6: iterations = "%sM" % (numIters / int(1e6))
	elif numIters >= 1e3: iterations = "%sK" % (numIters / int(1e3))
	else: iterations = "%s" % numIters

	prefix += ".k%s.%s.%s" % ("".join(map(str, ks)), iterations, acc)
	return prefix

def output_comet(args, mutations, results, collections, ks, runNum, maxPermutedWeight, permutedN):
	weight_func_mapping = {0: 'E', 1:'E', 2:'B', 3:'P'}
	header = "#Freq\tTotal Weight\tTarget Weight\t"
	header += "\t".join(["Gene set %s (k=%s)\tPhi %s\tWeight function %s" % (i, ks[i-1], i, i) for i in range(1, len(ks)+1)])
	tbl = [header]

	for S in collections:
		data = results[S]
		row = [ data["freq"], data["total_weight"], format(data["target_weight"], 'g') ]
		for d in sorted(data["sets"], key=lambda d: d["W"]):
			row += [", ".join(sorted(d["genes"])), d["prob"], weight_func_mapping[d["num_tbls"]] ]
		tbl.append("\t".join(map(str, row)))

	ensure_dir(args.output_dir)

	outputDirResults = args.output_dir + "/results/"
	outputDirViz = args.output_dir + "/viz/"
	ensure_dir(outputDirResults)
	ensure_dir(outputDirViz)

	# output results 
	fPrefix =  outputDirResults + iter_num('comet', runNum, ks, args.accelerator)	
	outputFile = "%s.tsv" % fPrefix	
	with open(outputFile, "w") as outfile: outfile.write( "\n".join(tbl) )
	paraJson = "%s.json" % fPrefix
	with open(paraJson, "w") as outfile: json.dump(vars(args), outfile)
	
	output_comet_viz(args, mutations, tbl, maxPermutedWeight, permutedN)

def output_comet_viz(args, mutations, resultsTable, maxPermutedWeight, permutedN):

	mutationMatrix   = args.mutation_matrix
	geneFile         = args.gene_file
	patientFile      = args.patient_file
	minFreq          = args.min_freq
	eventNamesFile   = args.event_names  	
	sampleToTypeFile = args.sample_types_file
	msf              = args.minimum_sampling_frequency 
	sec              = args.standard_error_cutoff 
	mew              = args.minimum_edge_weight 
	vizOutput        = args.output_dir +  "/viz"

	# Load the mutation data
	#mutations  = C.load_mutation_data(mutationMatrix, patientFile, geneFile, minFreq)
	m, n, genes, patients, geneToCases, patientToGenes = mutations
	eventNames = load_event_names(eventNamesFile, genes)

	###########################################################################
	# Construct marginal probability graph from the input CoMEt results file
	if args.verbose: print "* Constructing marginal probability graph..."

	res = construct_mp_graph( resultsTable, eventNames, msf, maxPermutedWeight )
	MPG, tables, expectedPoint = res
	edges = MPG.edges(data=True)

	if len(edges) == 0: # no significant results
		print "No significant collection. "
		exit(1)
		
	if args.verbose: print "\t- Edges:", len(edges)

	# Choose delta (the minimum edge weight in the marginal probability 
	# graph ) using a heuristic approach that selects delta at first elbow
	# with slope change > 0 using linear regression
	if args.verbose: print "* Choosing delta..."

	deltas = sorted(set( d['weight'] for u, v, d in MPG.edges(data=True)))
	realEdgeDist = compute_edge_dist(MPG, deltas)
	deltaPoint, edgeno = choose_delta(deltas, realEdgeDist, expectedPoint, sec)

	if args.verbose: print "\t- Delta: ", deltaPoint

	###########################################################################
	# Create the web output
	plot=None

	# Write the delta plot to file as an SVG, then load
	# it so we can embed it in the web page
	if args.verbose: print "* Plotting delta curve..."
	tmp = tempfile.mktemp(".svg", dir=".", prefix=".tmp")
	delta_plot(permutedN, deltas, realEdgeDist, tmp, expectedPoint, deltaPoint, edgeno)

	# Read in the graph, skipping the first four lines that are
	# extraneous header information that will only confuse a webpage
	with open(tmp) as f: plot = "".join(f.readlines()[4:])
	os.unlink(tmp)
	stats = dict(deltas=[dict(delta=deltaPoint, pval=0.)], plot=plot, N=permutedN)

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
	htmlOutput = "{}/index.html".format(vizOutput)
	with open(args.template_file) as template, open(htmlOutput, "w") as outfile:
		jsonData = json.dumps( dict(graph=graphData, mutations=mutations, tables=tables, stats=stats))
		html = template.read()
		html += "\n<script>\nvar data = {};\ndplusViz(data);\n</script>\n".format(jsonData)
		outfile.write( html )

	# Then copy the required JS and CSS files
	import shutil
	shutil.copyfile("comet/src/js/comet-viz.js", "{}/comet-viz.js".format(vizOutput))
	shutil.copyfile("comet/src/js/mp-graph.js", "{}/mp-graph.js".format(vizOutput))
	shutil.copyfile("comet/src/js/bower.json", "{}/bower.json".format(vizOutput))
	shutil.copyfile("comet/src/css/style.css", "{}/style.css".format(vizOutput))

###############################################################################
# Delta selection functions

def delta_plot(N, deltas, realEdgeDist, outfile, boundEdges, deltaPoint, edgeno):
	import matplotlib
	matplotlib.use('Agg') # for users without DISPLAY environment variable
	import matplotlib.pyplot as plt

	plt.rc('text', usetex=True)
	c1, c2, c3, c4 = [(0.9677975592919913, 0.44127456009157356, 0.5358103155058701), (0.5920891529639701, 0.6418467016378244, 0.1935069134991043), (0.21044753832183283, 0.6773105080456748, 0.6433941168468681), (0.6423044349219739, 0.5497680051256467, 0.9582651433656727)]
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
	""" Create the distribution of edge weights """
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

def choose_delta( deltas, realDist, expectedPoint, stdCutoff):
	""" Find the delta through L-elbow method """
	logY = np.log10(np.array(realDist))
	logX = np.log10(np.array(deltas))

	indices = -1
	maxDiff = -1

	# reverse order from large edge weight, set start_i as index passing expected edges
	startI = 0.
	for i in range(len(logX)-1, 0, -1):
		if realDist[i] > expectedPoint:
			startI = i
			break

	if len(logX) <= 3: # less and equal than three edge weights => can't do regression method. Output smallest edge.
		print "At most three edge weights in the marginal probability graph. Using the smallest edge weight as Delta..."
		return deltas[0], realDist[0]

	lastSlope = 0.
	i = startI # i: change point
	lastSlopeI = None # for storing the index of elbow point
	while True:
		iS = 0
        # define a line with standard mean error <= stdCutoff
		for j in range(2, len(logX)): 
			iS = i-j+1
			iE = i
			# linear regression
			regression = np.polyfit(logX[iS:iE+1] , logY[iS:iE+1], 1)
			# y=ax+b a:point_on_x, b: point_on_y
			pointOnX = regression[1] / (0 - regression[0])
			pointOnY = regression[1]
			# form line 
			line = regression[1] + regression[0] * logX[iS:iE+1]
			lineX = (logY[iS:iE+1]-regression[1])/(regression[0])
			# calculate standard mean error
			errX=math.sqrt(sum((lineX-logX[iS:iE+1])**2)/(iE+1-iS))
			errY=math.sqrt(sum((line-logY[iS:iE+1])**2)/(iE+1-iS))

			# if iS goest to the smallest edge weight or err of line > stdCutoff
			if (errX + errY)/2 > stdCutoff or iS == 0: 
				break

		if iS == 0 and lastSlopeI == None: # no elbow point
			print "No elbow point found! Try to lower the standard error cutoff with -rmse < ", stdCutoff
			exit(1)

		
		if lastSlopeI == None: # form first line. Assign slope and index for comparing with another line in the next round
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

def nCr(n,r):
	""" n choose r """
	return factorial(n) / factorial(r) / factorial(n-r)

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

###############################################################################
# Construct marginal prob. graph

def construct_mp_graph(resultsTable, eventNames, minSamplingFreq, maxPWeight):
	""" Construct marginal prob. graph from collections with weight larger than permuted data (maxPWeight) """
	""" Using genesets with the maximum weight to generate expected number of nodes in the graph for delta selection """
	tables = []
	geneSets, freqs = [], []
	N = 0
	
	maxRWeight = 0
	expectedPoint = 0
	
	for l in resultsTable:
		if l.startswith("#") or float(l.split("\t")[1]) < maxPWeight: continue
		if float(l.split("\t")[1]) > maxRWeight:
			maxRWeight = float(l.split("\t")[1])
			expectedPoint = sum([nCr(len(G.split(", ")),2) for G in l.rstrip().split("\t")[3::3] ])

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

	return MPG, tables, expectedPoint


def extract_relevant_edges(edges, minEdgeWeight):
	# Sort edges ascending by weight
	edges = [ (u, v, d) for u, v, d in edges if d['weight'] >= minEdgeWeight ]

	# Make a dictionary mapping 
	weightToEdges = defaultdict(list)
	for u, v, d in edges:
		weightToEdges[d['weight']].append( (u, v) )

	return sorted(weightToEdges.keys(), reverse=True), weightToEdges



