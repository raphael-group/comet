#!/usr/bin/python

###############################################################################
# Functions for loading mutation data
def load_mutation_data(filename, patientFile=None, geneFile=None, minFreq=0):
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
    """
    # Load the whitelists (if applicable)
    if patientFile:
        with open(patientFile) as f:
            patients = set( l.rstrip().split()[0] for l in f if not l.startswith("#") )
    else:
        patients = None

    if geneFile:
        with open(geneFile) as f:
            genes = set( l.rstrip().split()[0] for l in f if not l.startswith("#") )
    else:
        genes = set()

    # Parse the mutation matrix
    from collections import defaultdict
    geneToCases, patientToGenes = defaultdict(set), defaultdict(set)
    with open(filename) as f:
        arrs = [ l.rstrip().split("\t") for l in f if not l.startswith("#") ]
        for arr in arrs:
            patient, mutations = arr[0], set(arr[1:])

            if not patients or patient in patients:
                if genes: mutations &= genes
                #else: genes |= mutations                

                patientToGenes[patient] = mutations
                for gene in mutations:
                    geneToCases[gene].add(patient)

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
    return m, n, genes, patients, geneToCases, patientToGenes


def adj_dict_to_lists(xs, ys, d):
    """Convert a dictionary of x -> y to a list of lists, where
       each x corresponds to a list of indices of y."""
    M = []    
    for x, y_list in d.iteritems():
        M.append( [ ys.index(y) for y in y_list ] )
    return M

def convert_mutations_to_C_format(m, n, genes, patients, geneToCases, patientToGenes, subtypes=None):
    """We convert the dictionaries to list of lists so they're easier to parse in C."""

    if subtypes: genes += subtypes
    geneToIndex = dict(zip(genes, range(m)))
    indexToGene = dict(zip(range(m), genes))
    iPatientToGenes = adj_dict_to_lists(patients, genes, patientToGenes)
    iGeneToCases = adj_dict_to_lists(genes, patients, geneToCases)
    geneToNumCases = [ len(geneToCases[g]) for g in genes ]
    return iPatientToGenes, iGeneToCases, geneToNumCases, geneToIndex, indexToGene
