#!/usr/bin/python

###############################################################################
# Functions for loading mutation data
def load_mutation_data(filename, patientFile=None, geneFile=None, minFreq=0, subtypeFile=None):
    """Loads the mutation data in the given file.

    :type filename: string
    :param filename: Location of mutation data file.
    :type patient_file: string
    :param patient_file: Location of patient (whitelist) file.
    :type gene_file: string
    :param gene_file: Location of gene (whitelist) file.
    :rtype: Tuple

    **Returns:**
      * **m** (*int*) - number of patients.
      * **n** (*int*) - number of genes.
      * **genes** (*list*) - genes in the mutation data.
      * **patients** (*list*) - patients in the mutation data.
      * **geneToCases** (*dictionary*) - mapping of genes to the patients in which they are mutated.
      * **patientToGenes** (*dictionary*) - mapping of patients to the genes they have mutated.
      * **subtypes** (*list*) - list of unique subtypes.
    """
    # Load the patient whitelist (if applicable)
    if patientFile:
        with open(patientFile) as f:
            patients = set( l.rstrip().split()[0] for l in f if not l.startswith("#") )
    else:
        patients = None

    # Load the subtype information (if applicable)
    from collections import defaultdict
    subtypeDict = defaultdict( lambda: None, dict() )
    subtypes = set()
    if subtypeFile:
        subtypeDict = defaultdict( lambda: None, dict() )
        subtypes = set()
        with open( subtypeFile ) as f:
            sts = [ l.rstrip().split( '\t' ) for l in f if not l.startswith( '#' ) ]
            for p, s in sts:
                subtypeDict[p] = s
                subtypes.add( s )

    # Load the gene whitelist (if applicable)
    if geneFile:
        with open(geneFile) as f:
            genes = set( l.rstrip().split()[0] for l in f if not l.startswith("#") )
        genes |= subtypes
    else:
        genes = set()

    # Parse the mutation matrix
    geneToCases, patientToGenes = defaultdict(set), defaultdict(set)
    with open(filename) as f:
        arrs = [ l.rstrip().split("\t") for l in f if not l.startswith("#") ]
        for arr in arrs:
            patient, mutations = arr[0], set(arr[1:])

            if subtypeDict[patient]:
                mutations = mutations.union( subtypes.difference( set( [subtypeDict[patient]] ) ) )

            if not patients or patient in patients:
                if genes: mutations &= genes
                #else: genes |= mutations

                patientToGenes[patient] = mutations
                for gene in mutations:
                    geneToCases[gene].add(patient)

    genes = geneToCases.keys()
    # Remove genes with fewer than min_freq mutations
    toRemove = [ g for g in genes if len(geneToCases[g]) < minFreq ]
    for g in toRemove:
        for p in geneToCases[g]:
            patientToGenes[p].remove(g)
        del geneToCases[g]
        genes.remove(g)

    # Format and return output
    genes, patients = geneToCases.keys(), patientToGenes.keys()
    m, n = len(genes), len(patients)
    return m, n, genes, patients, geneToCases, patientToGenes, subtypes


def adj_dict_to_lists(xs, ys, d):
    """Convert a dictionary of x -> y to a list of lists, where
       each x corresponds to a list of indices of y."""
    M = []
    for x, y_list in d.iteritems():
        M.append( [ ys.index(y) for y in y_list ] )
    return M

def convert_mutations_to_C_format(m, n, genes, patients, geneToCases, patientToGenes, subtypes=None):
    """We convert the dictionaries to list of lists so they're easier to parse in C."""

    if subtypes:
        newg = set(genes).difference(set(subtypes))
        genes = list(newg)
        for s in subtypes:
            genes.append(s)
        #genes += subtypes

    geneToIndex = dict(zip(genes, range(m)))
    indexToGene = dict(zip(range(m), genes))
    iPatientToGenes = adj_dict_to_lists(patients, genes, patientToGenes)
    iGeneToCases = adj_dict_to_lists(genes, patients, geneToCases)
    geneToNumCases = [ len(geneToCases[g]) for g in genes ]
    return iPatientToGenes, iGeneToCases, geneToNumCases, geneToIndex, indexToGene
