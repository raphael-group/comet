#!/usr/bin/env python

# Import items from local modules
from mutation_data import *
from convergence import *
from constants import *
from permute import *

# There is a scipy stats function called binom_test, so we need to import
# the CoMET C modules last
from cComet import *

# Permutation test, conditioning on the frequency of each mutation
def py_permutation_test(k, N, tbl, num_permutations):
    # Compute the marginal frequencies per gene
    M = np.zeros((k, N), dtype=bool)
    sample_index = tbl[0]
    for i in range(1, 2 ** k):
        indices = [ j for j in range(k) if (1 << j) & i != 0 ]
        M[indices, sample_index:sample_index+tbl[i]] = True
        sample_index += tbl[i]

    # Compute's the number of samples with an exclusive mutation
    def T(A): return np.sum(A.sum(axis=0) == 1)

    # Create the permuted datasets and compute the significance
    obs = T(M)
    more_extreme_count = 0.
    for _ in range(num_permutations):
        for i in range(k): np.random.shuffle(M[i])
        if T(M) >= obs: more_extreme_count += 1.

    return more_extreme_count/num_permutations

def phi(k, N, tbl, exact_pvalthresh=0.001, binom_pvalthresh=0.005,
        co_cutoff=10, num_permutations=1000):
    # We use the exact test for k <= 3
    if k <= 3:
        return exact_test(k, N, tbl, 1.1)[1], EXACT
    # Otherwise we run the pipeline
    else:
        # Compute the binomial approximation and use it if (a) there are a
        # minimum number of co-occurrences (since the approximation is poor
        # for mostly exclusive cases) and (b) the P-value is above a minimum
        # value
        binom_pval = binom_test(k, N, list(tbl), 1.1)
        conum = sum(tbl[1:]) - sum(tbl[2**i] for i in range(k))
        if conum > co_cutoff or binom_pval > binom_pvalthresh:
            return binom_pval, BINOM

        # If the binomial doesn't work, then try the exact test. The exact test
        # runs until either the P-value is computed, or accumulated mass
        # reaches a certain value.
        num_tbls, exact_pval = exact_test(k, N, tbl, 1.1)
        if exact_pval != -1: # the P-value was computed exactly
            return exact_pval, EXACT
        else: # the P-value never finished
            return py_permutation_test(k, N, tbl, num_permutations), PERMUTATIONAL
