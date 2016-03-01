
###############################################################################
# Convert sampling solns to the output format
def convert_solns(indexToGene, solns):
    newSolns = []
    for arr in solns:
        arr.sort(key=lambda M: M[-2], reverse=True)
        S = tuple( frozenset([indexToGene[g] for g in M[:-2] ]) for M in arr )
        W = [ M[-2] for M in arr ]
        F = [ M[-1] for M in arr ]
        newSolns.append( (S, W, F) )

    return newSolns

###############################################################################
# Generate initial solutions 

# Call multi-dendrix to obtaining initial solns
def call_multidendrix(multi_dendrix, mutations, k, t):
    alpha, delta, lmbda = 1.0, 0, 1 # default of multidendrix
    geneSetsWithWeights = multi_dendrix.ILP( mutations, t, k, k, alpha, delta, lmbda)    
    multiset = list()
    for geneSet, W in geneSetsWithWeights:
        for g in geneSet:
            multiset.append(g)
    return multiset

# Generate solns based on user-input.
def initial_solns_generator(r, mutations, ks, assignedInitSoln, subtype, importMultidendrix, multi_dendrix):
    runInit = list()
    totalOut = list()

    # Assign input solns
    if assignedInitSoln:
        if len(assignedInitSoln) == sum(ks):
            #print 'load init soln', "\t".join(assignedInitSoln)
            runInit.append(assignedInitSoln)
        elif len(assignedInitSoln) < sum(ks): # fewer initials than sampling size, randomly pick from genes
            import random
            rand = assignedInitSoln + random.sample(set(mutations[2])-set(assignedInitSoln), sum(ks)-len(assignedInitSoln))
            #print 'load init soln with random', rand
        else:
            sys.stderr.write('Too many initial solns for CoMEt.\n')
            exit(1)

    # Load Multi-Dendrix results as initialization 
    if importMultidendrix and not subtype and ks.count(ks[0])==len(ks):

        md_init = call_multidendrix(multi_dendrix, mutations, ks[0], len(ks))
        #print ' load multi-dendrix solns', md_init
        runInit.append(list(md_init))

    # assign empty list to runInit as random initials
    for i in range(len(runInit), r):
        runInit.append(list())

    for i in range(r):
        totalOut.append(dict())

    return runInit, totalOut

# Load precomputed scores from input file  (--precomputed_scores)
def load_precomputed_scores(infile, mutations, subt):

    if subt: mutations = mutations + (subt,)
    cMutations = C.convert_mutations_to_C_format(*mutations)
    iPatientToGenes, iGeneToCases, geneToNumCases, geneToIndex, indexToGene = cMutations

    baseI = 3  # sampling freq., total weight, target weight
    setI = 3 # gene set, score, weight function

    matchObj = re.match( r'.+\.k(\d+)\..+?', infile)

    loadingT = len(matchObj.group(1)) # determine t:the number of gene sets.
    for l in open(infile):
        if not l.startswith("#"):
            v = l.rstrip().split("\t")
            j = 0
            for i in range(loadingT):
                gSet = [geneToIndex[g] for g in v[baseI + j].split(", ")]
                C.load_precomputed_scores(float(v[baseI + j + 1]), len(v[baseI + j].split(", ")), int(v[baseI + j + 2]), gSet)
                j += setI

###############################################################################
# Merge results
def merge_results(convResults):
    total = dict()
    for results in convResults:
        for key in results.keys():
            if key in total:
                total[key]["freq"] += results[key]["freq"]
            else:
                total[key] = results[key]

    return total

def merge_runs(resultsPre, resultsNew):

    for key in resultsNew.keys():
        if key in resultsPre:
            resultsPre[key]["freq"] += resultsNew[key]["freq"]
        else:
            resultsPre[key] = resultsNew[key]
