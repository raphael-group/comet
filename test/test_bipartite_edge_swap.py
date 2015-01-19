#!/usr/bin/python

import networkx as nx, sys, os, time, random, numpy as np
sys.path.append("../comet")
from permute_matrix import bipartite_edge_swap
from networkx.algorithms import bipartite
import permute as P

def check_degrees(G, H):
    if len(set(G.nodes()) ^ set(H.nodes())) != 0:
        print set(G.nodes()) - set(H.nodes())
        print set(H.nodes()) - set(G.nodes())
    for i, n in enumerate(G.nodes()):
        if H.degree(n) != G.degree(n):
            print i, n, H.degree(n), G.degree(n)
            return False
    return True

def call_bipartite_edge_swap( G, xs, ys, Q ):
    # Compute the desired pieces of the graph structure
    degrees = [ G.degree(n) for n in xs + ys ]
    A = np.array(bipartite.biadjacency_matrix(G, row_order=xs, column_order=ys, dtype=np.int32))
    
    # Set up and call the permute matrix function
    max_tries = 1e75
    seed      = time.time()
    nswap     = len(G.edges()) * Q
    B = bipartite_edge_swap(A, degrees, nswap, max_tries)
    H = nx.Graph()
    H.add_edges_from([ (xs[u], ys[v]) for u, v in zip(*np.where(B == 1)) ])
    return H

def test_random_bipartite_graph(numXs, numYs, numEdges ):
    # Create a random bipartite graph
    G = nx.Graph()
    genes = ["G{}".format(i) for i in range(1, numXs+1) ]
    patients = ["TCGA-{}".format(i) for i in range(1, numYs+1) ]
    assert(numEdges < numXs * numYs)
    xs, ys = set(), set()
    for i in range(numEdges):
        x = random.choice(genes)
        y = random.choice(patients)
        G.add_edge(x,y)
        xs.add(x)
        ys.add(y)
    return list(xs), list(ys), G

if __name__ == "__main__":
    Q = 100
    for i, (numXs, numYs, numEdges) in enumerate([(300, 400, 20000)]):
        # Create a test graph
        xs, ys, G = test_random_bipartite_graph(numXs, numYs, numEdges)
        print "Test {}: {} genes x {} samples, {} mutations".format(i+1, len(xs), len(ys), numEdges)
        
        # Run and time in Fortran
        start = time.time()
        H = call_bipartite_edge_swap( G, xs, ys, Q=Q )
        print "\tFortran:", time.time() - start, 'secs'
        worked = check_degrees(G, H)
        if not worked:
            raise ValueError("Degrees in permuted graph are different than original graph.")
        
        # Run and time in Python
        start = time.time()
        H = P.bipartite_double_edge_swap( G, xs, ys, nswap=numEdges * Q)
        print "\tPython:", time.time() - start, 'secs'
        check_degrees(G, H)
        if not worked:
            raise ValueError("Degrees in permuted graph are different than original graph.")

        print

