#include "cometmodule.h"

/* Keeps track of all fields during MCMC */
/*******************************************************************************
  Contingency Table Generation
 *******************************************************************************/
void contingency_tbl_gen(mutation_data_t *mut_data, int k, int *ctbl, int *gene_set) {
    int buf;
    int index, layer, ones;
    int i, j; 
    int tmp[k];
    int ctbl_len = 1 << k;

    for (i=0; i<k; i++) {
        tmp[i] = gene_set[(k-1)-i];
    }

    // Initialize the contingency table to zeros
    memset(ctbl, 0, sizeof(int) * ctbl_len);   

    ones = (1 << k) - 1;
    for (i=0; i<mut_data->numPatients; i++) {
        buf = 0 << k;
        for (j=0; j < k; j++) {
            index = mut_data->P2G[i][tmp[j]];
            layer = ((index == (char) 1) ? ones : 0) & (1 << j);
            buf = buf | layer;
        }
        ctbl[buf] += 1;
    }    
}


// Loads mutation data into a matrix (etc) from PyObjects
/* Initialize dendrix parameters.
   P2numMuts will be an array where the ith entry holds the number of genes mutated in patient i. 
   P2G will be a 2-D array where the jth entry of the ith row is a bit indicating whether or not gene[j]
   G2numMutated will be an array where the ith entry is the total number of mutations in gene[i] 
   is mutated in patient i. 
 */
void init_mutation_data(mutation_data_t *mut_data, int numPatients, int numGenes,
                        PyObject *P2mutatedGenes, PyObject *gene2numMutations){
    // Create log factorial lookup table
    int i, j, p2g_index;
    int num_genes;
    PyObject *genes, *temp;
    mut_data->numPatients = numPatients;
    mut_data->numGenes = numGenes;

    for (i=0; i<mut_data->numPatients; i++) {
        genes = PyList_GetItem(P2mutatedGenes, i);
        num_genes = PyList_Size(genes);
        mut_data->P2numMuts[i] = num_genes;
        mut_data->P2G[i] = malloc(sizeof(char) * mut_data->numGenes);
        memset(mut_data->P2G[i], 0, sizeof(char) * mut_data->numGenes);
        for (j=0; j<num_genes; j++) {
            temp = PyList_GetItem(genes, j);
            p2g_index = (int) PyLong_AsLong(temp);
            mut_data->P2G[i][p2g_index] = (char) 1;
        }
    }
    for (i=0; i < mut_data->numGenes; i++) {
        mut_data->G2numMuts[i] = (int) PyLong_AsLong (PyList_GetItem(gene2numMutations, i));
    }
}

/* Allocate all memory needed for a dendrix parameter object */
mutation_data_t *mut_data_allocate(int numPatients, int numGenes) {
    mutation_data_t *mut_data  = (mutation_data_t*) malloc(sizeof(mutation_data_t));
    mut_data->P2G             = malloc(sizeof(char *) * numPatients);
    mut_data->P2numMuts       = malloc(sizeof(int) * numPatients);
    mut_data->G2numMuts       = malloc(sizeof(int) * numGenes);
    
    return mut_data;
}

/* Allocate all memory needed for a dendrix parameter object in permutation test*/
mutation_data_t *permute_mut_data_allocate(int numPatients, int numGenes) {
    mutation_data_t *mut_data  = (mutation_data_t*) malloc(sizeof(mutation_data_t));
    mut_data->P2G             = malloc(sizeof(char *) * numPatients);
    mut_data->P2numMuts       = calloc(numPatients, sizeof(int));
    mut_data->G2numMuts       = calloc(numGenes, sizeof(int));
    
    return mut_data;
}

/* Free all memory associated with a dendrix parameter object.*/
void mut_data_free(mutation_data_t *mut_data) {
    free_ptr_array((void **) mut_data->P2G, mut_data->numPatients);
    free(mut_data->P2numMuts);
    free(mut_data->G2numMuts);
    free(mut_data);
}

