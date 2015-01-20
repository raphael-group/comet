#!/usr/bin/python

# Import globally required modules
import sys, random, numpy as np
try:
        import networkx as nx
        from networkx.algorithms.bipartite import biadjacency_matrix
except ImportError:
	print 'Error!'
	print '\tCould not import NetworkX (http://networkx.github.com).'
	print '\tMake sure NetworkX is in your path.'
	sys.exit(1)

# Try to import the Fortran double edge swap, otherwise
# import the Python version
try: 
        from permute_matrix import bipartite_edge_swap
        fortranBindings = True
except ImportError:
        fortranBindings = False
        sys.stderr.write("[Warning] Could not import Fortran bipartite_edge_swap bindings.\n")

def bipartite_double_edge_swap(G, genes, patients, nswap=1, max_tries=1e75):
	"""A modified version of the double_edge_swap function in NetworkX to preserve the bipartite structure of the graph.
	
	For more details on this function, please see the `original NetworkX function <http://goo.gl/wWxBD>`_ that I shamelessly used to make this one.
	The only major change here is that I ensure that u,v and x,y are of the same type (i.e. genes or patients).
	"""
	if nswap>max_tries:
		raise nx.NetworkXError("Number of swaps > number of tries allowed.")
	if len(G) < 4:
		raise nx.NetworkXError("Graph has less than four nodes.")
	# Instead of choosing uniformly at random from a generated edge list,
	# this algorithm chooses nonuniformly from the set of nodes with
	# probability weighted by degree.
	n=0
	swapcount=0
	keys,degrees=zip(*G.degree().items()) # keys, degree
	cdf=nx.utils.cumulative_distribution(degrees)  # cdf of degree
	while swapcount < nswap:
		# pick two random edges without creating edge list
		# choose source node indices from discrete distribution
		(ui,xi)=nx.utils.discrete_sequence(2,cdistribution=cdf)
		if ui==xi:
			continue # same source, skip
		u=keys[ui] # convert index to label
		x=keys[xi]
		if (u in genes and x in genes) or (u in patients and x in patients):
			continue # both are genes, skip

		patient1 = u if u in patients else x
		gene1    = x if x in genes else u
		
		# choose target uniformly from neighbors
		patient2=random.choice( list(G[gene1]) )
		gene2=random.choice( list(G[patient1]) )
	
		# don't create parallel edges
		if (gene1 not in G[patient1]) and (gene2 not in G[patient2]): 
			G.add_edge(gene1,patient1)
			G.add_edge(gene2,patient2)
			
			G.remove_edge(gene1,patient2)
			G.remove_edge(patient1, gene2)
			swapcount+=1

		if n >= max_tries:
			e=('Maximum number of swap attempts (%s) exceeded '%n +
			   'before desired swaps achieved (%s).'%nswap)
			raise nx.NetworkXAlgorithmError(e)
		n+=1
	return G

def construct_mutation_graph(geneToCases, patientToGenes):
	genes, patients = geneToCases.keys(), patientToGenes.keys()
	nodes = genes + patients
	edges = [ (gene, patient) for gene in genes for patient in geneToCases[gene] ]
	G = nx.Graph()
	G.add_nodes_from(nodes)
	G.add_edges_from(edges)
	return G

def graph_to_mutation_data(H, genes, patients):
	geneToCases, patientToGenes = dict([(g, set()) for g in genes]), dict( )
	for patient in patients:
		mutations = H[patient]
		patientToGenes[patient] = set( mutations )
		for g in mutations: geneToCases[g].add( patient )
				
	genes, patients = geneToCases.keys(), patientToGenes.keys()
	m, n = len(genes), len(patients)
	return m, n, genes, patients, geneToCases, patientToGenes

def permute_mutation_data(G, genes, patients, seed, Q=100):
        if fortranBindings:
                # Compute the desired pieces of the graph structure
                xs = sorted(genes, key=G.degree, reverse=True)
                ys = sorted(patients, key=G.degree, reverse=True)
                x_degrees = [ G.degree(x) for x in xs ]
                y_degrees = [ G.degree(y) for y in ys ]
                A = np.array(biadjacency_matrix(G, row_order=xs,
                                                column_order=ys,
                                                dtype=np.int32))
                
                # Set up and call the permute matrix function
                B = bipartite_edge_swap(A, x_degrees, y_degrees, len(G.edges()) * Q, 1e75, seed=seed)
                H = nx.Graph()
                H.add_edges_from([ (xs[u], ys[v]) for u, v in zip(*np.where(B == 1)) ])
        else:
                H = G.copy()
                random.seed(seed)
                bipartite_double_edge_swap(H, genes, patients, nswap=Q * len( G.edges() ))
	return graph_to_mutation_data(H, genes, patients)
