#ifndef COMETMODULE_H
#define COMETMODULE_H

#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <glib.h>
#include <string.h>
#include <stdio.h>

#include "utils/uthash.h"

/* Keeps track of all fields during MCMC */
typedef struct mutation_data {
//    double (*weight_func) (struct mutation_data *self);
    int numGenes, numPatients;
    int *P2numMuts, *G2numMuts;
    char **P2G;
} mutation_data_t;

typedef struct current_values {
    int *current_soln, *current_table;
    double current_prob, current_weight;
    int current_num_tbls;
} current_soln_values_t;

/* Snapshots of all fields during MCMC */
typedef struct frozen_arrays {
    double *frozenadjweights;
    int **frozensets;
    double *frozenprobs;
    int **frozentables;
} frozen_arrays_t;

/* Contains fields associated with one gene set. Used in memoization.*/
typedef struct memo_data {
    int k;
    double prob;
    double weight;
    int num_tbls;
    int *table;
} memo_data_t;

// Structs for the hash
typedef struct{
    int genes[10];
} geneset_t;

typedef struct {
    geneset_t id; /* REQUIRED: the key */
    double weight;
    int function;
    UT_hash_handle hh; /* REQUIRED: makes this structure hashable */
} weight_t;

// Global store of the weights of seen gene sets
extern weight_t *weightHash;

/* The weight function used by dendrix. */
typedef double (*weight_function) (mutation_data_t *mut_data, current_soln_values_t *curr_soln_vals, int *ctbl, int k, double pvalthresh); 
typedef double (*target_function) (double weight, double c); 

extern char weight_func_c;
//extern mutation_data_t *glbl_mut_data;
//extern int log_results;

/*******************************************************************************************************/

/***************/
/* utilities.c */
/***************/
extern double *lnfacs; 
double lnfac(int a);
int min(int x, int y);
int max(int x, int y);
int arr_min(int *arr, int len);
int arr_max(int *arr, int len);
void precompute_factorials(int N);
void int_arr_to_string(int *array, int len, char *dest, char *delimiter); 
void convert_int_array_to_pylist(int len, PyObject *pyArr, int *arr);
void convert_double_array_to_pylist(int len, PyObject *pyArr, double *arr);
int indexOf(int n, int *array, int k);

//extern int *glbl_G2numMuts; 
/* Global array of genes to number of mutations for sorting */
int comp_g2s_numMutations(const void *e1, const void *e2);

int sum_int_array(int *arr, int len); 

void freeze_int_array(int i, int len, int **frozenarrays, int *arr); 
current_soln_values_t *current_arrays_allocate(int k); 
void current_arrays_free(current_soln_values_t *curr_soln_vals); 

void key_destroyed(gpointer gstr_key); 
void value_destroyed(gpointer data_value); 

frozen_arrays_t *frozen_allocate(int numFrozen); 
void frozen_free(int numFrozen, frozen_arrays_t *frozen); 

void free_ptr_array(void **array, int len);
void frozensets2pylists(frozen_arrays_t *frozen, int numFrozen, int k, PyObject *solns, PyObject *weights, PyObject *tables, PyObject *probs); 

void copyArr(double *A, double *B, int n);
double sum(double *arr, int n);
int sum_array(int a[], int start, int end);
int ascending(const void * elem1, const void * elem2);

/*********************/
/* weights.c */
/*********************/
double exact_test(int k, int numSamples, int *ctbl, int *num_tbls);

double comet_exact_test(int k, int N, int *ctbl, int *final_num_tbls, double pvalthresh);
double comet_binomial_test(int k, int N, int *tbl, double pvalthresh);
double comet_binomial_co_test(int k, int N, int *tbl, double pvalthresh);
double comet_permutation_test(int k, int num_permutation, int num_patients, int *gene_set, int *freq, int *tbl);
double binomial_cdf(int obs, int num_trials, double p);
void comet_phi(int *genes, int k, int n, mutation_data_t *A, int *ctbl, int** kelem, double *score, int *func
                , int co_cutoff, int permutation_iteration
                , double binom_pvalthreshold, double pvalthresh);

double dendrix(int *ctbl, int k);

int *get_ex_cells(int k);
int *get_co_cells(int k);
int num_ones(int n);
int sum_cells( int *tbl, int *cells, int num_entries);

/*******************/
/* mutation_data.c */
/*******************/
void contingency_tbl_gen(mutation_data_t *mut_data, int k, int *ctbl, int *soln); 
void init_mutation_data(mutation_data_t *mut_data, int numPatients, int numGenes,
                        PyObject *P2mutatedGenes, PyObject *gene2numMutations);
/* Allocate all memory needed for a dendrix parameter object */
mutation_data_t *mut_data_allocate(int numPatients, int numGenes); 
mutation_data_t *permute_mut_data_allocate(int numPatients, int numGenes); 
/* Free all memory associated with a dendrix parameter object.*/
void mut_data_free(mutation_data_t *mut_data); 


/**********************/
/* comet_exhaustive.c   */
/**********************/

frozen_arrays_t *comet_exhaustive(mutation_data_t *mut_data, int k, int *numSets, double pvalthresh); 
int choose(int N, int k); 

/* Functions callable from Python*/
PyObject *py_free_factorials(PyObject *self, PyObject *args);
PyObject *py_exact_test(PyObject *self, PyObject *args);
PyObject *py_binomial_test(PyObject *self, PyObject *args);
PyObject *py_comet(PyObject *self, PyObject *args);
PyObject *py_comet_score(PyObject *self, PyObject *args);
PyObject *py_exhaustive(PyObject *self, PyObject *args);
PyObject *py_load_precomputed_scores(PyObject *self, PyObject *args);

/**********************/
/* comet_mcmc.c   */
/**********************/

void comet_mcmc(mutation_data_t *A, int m, int n, int *ks, int t, int num_iters, int step_len, int amp, 
                    int nt, double binom_cutoff, int *initial_soln, int size_initial_soln,
                    int size_subtype, double pval_cutoff, int **gene_sets, double **set_weights, int **functions, 
                    int verbose);

#endif
