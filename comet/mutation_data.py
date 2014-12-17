###############################################################################
# Functions for loading mutation data
def load_mutation_data_w_cutoff(file_loc, patient_wlst=None,
                                gene_wlist=None, cutoff=0):
    """Loads the mutation data in the given file, restricting to genes with a given mutation frequency. 

    :type file_loc: string
    :param file_loc: Location of mutation data file.
    :type patient_wlst: dictionary
    :param patient_wlst: Maps patient IDs to whether they should be included in the analyzed mutation data.
    :type gene_wlist: dictionary
    :param gene_wlist: Maps genes to whether they should be included in the analyzed mutation data.
    :type cutoff: int
    :param cutoff: Minimum mutation frequency a gene must have to be included in the analyzed mutation data.
    :returns: Mutation data tuple (see :func:`load_mutation_data`).

    **Example:**
      A view into the example data:
        >>> file_loc = 'test.m2'
        >>> open(file_loc).readlines() # view of the  data
        ["TCGA-01\\tG1\\tG3\\tG5\\n", "TCGA-02\\tG2\\tG1\\tG4\\n"]
    
      Mutation data with no cutoff:
        >>> load_mutation_data_w_cutoff(file_loc, cutoff=0) # mutation data with no cutoff
        (5, 2, ["G1", "G2", "G3", "G4", "G5"], ["TCGA-01", "TCGA-02"],
            {"G1" : ["TCGA-01", "TCGA-02"], "G2" : ["TCGA-02"], "G3" : ["TCGA-01"],
            "G4" : ["TCGA-2"], "G5" : ["TCGA-01"]},
            {"TCGA-01" : ["G1", "G3", "G5"], "TCGA-02" : ["G2", "G1", "G4"]})
      
      Mutation data with a cutoff of 2:
        >>> load_mutation_data_w_cutoff(file_loc, cutoff=2)
        (1, 2, ["G1"], ["TCGA-01", "TCGA-02"], {"G1" : ["TCGA-01", "TCGA-02"]},
            {"TCGA-01" : ["G1"], "TCGA-02" : ["G1"]})
    
    **See also:**
    :func:`white_and_blacklisting`, :func:`load_mutation_data`

    """
    mutation_data = load_mutation_data(file_loc, patient_wlst, gene_wlist)
    if cutoff == 0: return mutation_data
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    for gene in genes:
        if len(mutation2patients[gene]) < cutoff:
            for patient in mutation2patients[gene]:
                patient2mutations[patient].remove( gene )
            del mutation2patients[gene]
    
    m, genes = len(mutation2patients.keys()), mutation2patients.keys()
    return m, n, genes, patients, mutation2patients, patient2mutations

def load_mutation_data(file_loc, patient_whitelist=None, gene_whitelist=None):
    """Loads the mutation data in the given file. 

    :type file_loc: string
    :param file_loc: Location of mutation data file.
    :type patient_whitelist: dictionary
    :param patient_whitelist: Maps patient IDs to whether they should be included in the analyzed mutation data.
    :type gene_whitelist: dictionary
    :param gene_whitelist: Maps genes to whether they should be included in the analyzed mutation data.
    :rtype: Tuple

    **Returns:**
      * **m** (*int*) - number of patients.
      * **n** (*int*) - number of genes.
      * **genes** (*list*) - genes in the mutation data.
      * **patients** (*list*) - patients in the mutation data.
      * **mutation2patients** (*dictionary*) - mapping of genes to the patients in which they are mutated.
      * **patient2mutations** (*dictionary*) - mapping of patients to the genes they have mutated.

    **Example:**
      A view into example data:
        >>> file_loc = 'test.m2'
        >>> open(file_loc).readlines()
        ["TCGA-01\\tG1\\tG3\\tG5\\n", "TCGA-02\\tG2\\tG1\\tG4\\n"]
      
      Mutation data with no whitelisting:
        >>> load_mutation_data(file_loc)
        (5, 2, ["G1", "G2", "G3", "G4", "G5"], ["TCGA-01", "TCGA-02"],
            {"G1" : ["TCGA-01", "TCGA-02"], "G2" : ["TCGA-02"], "G3" : ["TCGA-01"],
            "G4" : ["TCGA-2"], "G5" : ["TCGA-01"]},
            {"TCGA-01" : ["G1", "G3", "G5"], "TCGA-02" : ["G2", "G1", "G4"]})
      
      Mutation data with patient whitelisting only:
        >>> patient_wlst = {"TCGA-01" : False, "TCGA-02" : True}
        >>> load_mutation_data_w_cutoff(file_loc, patient_wlst)
        (3 , 1, ["G1", "G2", "G4"], ["TCGA-02"],
            {"G1" : ["TCGA-02"], "G2" : ["TCGA-02"], "G4" : ["TCGA-01"]},
            {"TCGA-01" : ["G1", "G2", "G4"]})
      
      Mutation data with patient and gene whitelisting:
        >>> gene_wlst = {"G1" : True, "G2" : True, "G4" : False}
        >>> load_mutation_data_w_cutoff(file_loc, patient_wlst, gene_wlst
        (2 , 1, ["G1", "G2"], ["TCGA-02"], {"G1" : ["TCGA-02"], "G2" : ["TCGA-01"]},
            {"TCGA-01" : ["G1", "G2"]})

    **See also:**
    :func:`white_and_blacklisting`, :func:`load_mutation_data_w_cutoff`
    """
    patient2mutations, mutation2patients = {}, {}
    for line in [l.rstrip() for l in open(file_loc) if not l.startswith('#')]:
        arr = filter(lambda g: g != "", line.split('\t'))
        patient, mutations = arr[0], arr[1:]

        if gene_whitelist:
            mutations = filter(lambda g: gene_whitelist[g], mutations)
        if patient_whitelist and not patient_whitelist[patient]: continue
        
        patient2mutations[patient] = mutations
        for gene in mutations:
            if not gene in mutation2patients.keys():
                mutation2patients[gene] = set([patient])
            else:
                mutation2patients[gene].add(patient)

    genes, patients = mutation2patients.keys(), patient2mutations.keys()
    m, n = len(genes), len(patients)
    return m, n, genes, patients, mutation2patients, patient2mutations

def white_and_blacklisting(patient_wlst=None, patient_blst=None, gene_wlst=None,
                           gene_blst=None):
    '''Reconciles the different white- and blacklists provided as input into Multi-Dendrix.

    :type patient_wlst: string
    :type patient_blst: string
    :type gene_wlst: string
    :type gene_blst: string
    :param patient_wlst: File location of patients to be *included* in analyzed mutation data.
    :param patient_blst: File location of patients to be *excluded* in analyzed mutation data.
    :param gene_wlst: File location of patients to be *included* in analyzed mutation data.
    :param gene_blst: File location of patients to be *excluded* in analyzed mutation data.

    **Returns:**
      * gene2include (*dictionary*): mapping of genes to whether they should be included in the analyzed mutation data.
      * patient2include (*dictionary*): mapping of patients to whether they should be included in the analyzed mutation data.

    **Examples:**
      *(For brevity, examples are for patient white- and blacklists only)*

      Patient whitelisting only:
        >>> patient_wlst = 'patient.wlst'
        >>> open(patient_wlst).readlines()
        ["TCGA-01", "TCGA-02", "TCGA-03"]
        >>> white_and_blacklisting(patient_wlst)
        (defaultdict(<function <lambda>>, {}), {"TCGA-01", "TCGA-02", "TCGA-03"})
      Conflicting patient white- and blacklists (whitelist has final word):
        >>> patient_blst = 'patient.blst'
        >>> open(patient_wlst).readlines()
        ["TCGA-02", "TCGA-04"]
        >>> white_and_blacklisting(patient_wlst)
        (defaultdict(<function <lambda>>, {}), {"TCGA-01", "TCGA-02", "TCGA-03"})


    **See also:** :func:`load_mutation_data`, :func:`load_mutation_data_w_cutoff`.

    '''

    # Blacklisting and whitelisting works as follows. If a whitelist is passed in,
    # only genes/patients in the whitelist and NOT in the blacklist are allowed. If
    # only a blacklist is passed in, all genes/patients not on the blacklist are
    # allowed.
    from collections import defaultdict
    if patient_wlst:
        patient2include = defaultdict(lambda : False)
        patient_whitelist = [l.rstrip().split()[0] for l in open(patient_wlst)
                             if not l.startswith("#")]

        for p in patient_whitelist: patient2include[p] = True
    else:
        patient_whitelist = None
        patient2include = defaultdict(lambda : True)

    if patient_blst:
        patient_blacklist = [l.rstrip() for l in open(patient_blst)]
        for p in patient_blacklist: patient2include[p] = False

    if gene_wlst:
        gene2include = defaultdict(lambda : False)
        gene_whitelist = set([l.split()[0] for l in open(gene_wlst)])
        for g in gene_whitelist: gene2include[g] = True
    else: gene_whitelist, gene2include = None, defaultdict(lambda : True)
    
    if gene_blst:
        gene_blacklist = [l.rstrip() for l in open(gene_blst)]
        for g in gene_blacklist: gene2include[g] = False

    return gene2include, patient2include

# Converts a mutation2patients hash into a 2D-array where the ith row
# contains the indices of patients that have mutations in genes[i]
def convert_mutation_data_to_int(patients, genes, mutation2patients):
    iG2P = []    
    for g in genes:
        mutations = []
        for p in mutation2patients[g]:
            mutations.append( patients.index(p) )
        iG2P.append( mutations )
    return iG2P

# Converts a patients2genes hash into a 2D-array where the ith row
# contains the indices of genes mutated in that patient.
def convert_p2m_to_int(patients, genes, patient2mutations):
    iP2G = []
    for patient in patient2mutations:
        mutated_genes = []
        for i in range(len(genes)):
            if genes[i] in patient2mutations[patient]:
                mutated_genes.append(i)
        iP2G.append(mutated_genes)
    return iP2G

def convert_mutation_data_for_c(genes, patients, mutation2patients, patient2mutations):
    # m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    # genes.sort(key=lambda g: len(mutation2patients[g]), reverse=True)
    gene2index = dict([ (genes[i], i) for i in range(len(genes)) ])
    index2gene = dict([ (i, genes[i]) for i in range(len(genes)) ])

    iP2G = convert_p2m_to_int(patients, genes, patient2mutations)
    iG2P = convert_mutation_data_to_int(patients, genes, mutation2patients)
    iG2numMuts = [len(mut_pat_list) for mut_pat_list in iG2P]
    return iP2G, iG2P, iG2numMuts, gene2index, index2gene

def convert_mutation_data_for_c_with_subtype(genes, patients, mutation2patients, patient2mutations, sub):
    # m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    # genes.sort(key=lambda g: len(mutation2patients[g]), reverse=True)    
    newg = set(genes).difference(set(sub))
    genes = list(newg)
    for s in sub:
        genes.append(s)

    gene2index = dict([ (genes[i], i) for i in range(len(genes)) ])
    index2gene = dict([ (i, genes[i]) for i in range(len(genes)) ])

    iP2G = convert_p2m_to_int(patients, genes, patient2mutations)
    iG2P = convert_mutation_data_to_int(patients, genes, mutation2patients)
    iG2numMuts = [len(mut_pat_list) for mut_pat_list in iG2P]
    return iP2G, iG2P, iG2numMuts, gene2index, index2gene
