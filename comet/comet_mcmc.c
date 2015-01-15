#include <Python.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>

#include "cometmodule.h"

/*****************************************************************************/
/* Computing the weight */
int lookups = 0;
int calculations = 0;
double pvalthresh = 0.001;
int mcmc_amp = 2;
int permutation_iteration = 1000;
int co_cutoff = 10;
double binom_pvalthreshold = 0.005;
bool dendrix_pass = false;

int compute_num_co(int *ctbl, int k){  
  int *co_cells, num_co_cells=pow(2, k)-k-1, numcells;
  co_cells = get_co_cells(k);  
  numcells = sum_cells(ctbl, co_cells, num_co_cells);
  free(co_cells);
  return numcells;
}

int getFunction(int *set, int t_i, int *ks, int *gi_sum){
  geneset_t *key;
  weight_t *w;
  int y, z;
  int start = gi_sum[t_i] - ks[t_i];
  int end = gi_sum[t_i];
  int k = end-start;
  // Generate the key by sorting the arr
  
  key = malloc( sizeof(geneset_t) );
  memset(key, 0, sizeof(geneset_t));
  for (y = start, z = 0; y < end; y++, z++) key->genes[z] = set[y];
  qsort(key->genes, k, sizeof(int), ascending);
  HASH_FIND(hh, weightHash, key, sizeof(geneset_t), w);  /* id already in the hash? */
  if (w == NULL){
    return 0;
  }
  else{
    return w->function;
  }

}
// Get the weight for a geneset, either from the hash
// or by computing it then adding it to the hash, start: start index of t-th gene set to exchange
double getWeight(int *set, int t_i, int *ks, int *gi_sum, int** kelem, int n, mutation_data_t *A){
  // Declarations
  geneset_t *key;
  weight_t *w;
  int i, j;
  int start = gi_sum[t_i] - ks[t_i];
  int end = gi_sum[t_i];
  int k = end-start;
  int num_entries = 1 << k;  
  int ctbl[num_entries];      
  
  // Generate the key by sorting the arr
  key = malloc( sizeof(geneset_t) );
  memset(key, 0, sizeof(geneset_t));
  for (i = start, j = 0; i < end; i++, j++) key->genes[j] = set[i];
  qsort(key->genes, k, sizeof(int), ascending);

  // Try and find the geneset
  HASH_FIND(hh, weightHash, key, sizeof(geneset_t), w);  /* id already in the hash? */

  // If it's not there, add it
  if (w == NULL) {
    calculations += 1;
    w = malloc(sizeof(weight_t));
    w->id = *key;
    //w->weight = W(start, k, n, A, set);
    contingency_tbl_gen(A, k, ctbl, key->genes); 

    if(dendrix(ctbl, k) > 0){
      double score;
      int func;
      comet_phi(key->genes, k, n, A, ctbl, kelem, &score, &func, co_cutoff, permutation_iteration, binom_pvalthreshold, pvalthresh);
      w->weight = -log(pow(score, mcmc_amp));
      w->function = func;
      HASH_ADD( hh, weightHash, id, sizeof(geneset_t), w );  /* id: name of key field */
      dendrix_pass = true;
    }
    else{      
      w->weight = 0.0;
      HASH_ADD( hh, weightHash, id, sizeof(geneset_t), w );  /* id: name of key field */
    }
  }
  else{
    lookups += 1;
    dendrix_pass = true;
  }

  // Return the gene set's weight
  return w->weight;
}

void delete_all(){
  weight_t *current, *tmp;
  HASH_ITER(hh, weightHash, current, tmp){
    HASH_DEL(weightHash, current);
    free(current);
  }

}

/*****************************************************************************/
// Generate a collection of t random sets of k genes (integers)
void rand_soln(int *set, double *weights, int m, int n, int t, int* gi_sum, int* ks, int** kelem, int num_genes, mutation_data_t *A){
  // Declarations
  int i, r, j=0;

  // Select t random gene sets of size k, making sure
  // there are no duplicates
  for(i = 0; i < num_genes; i++){
    do{
      r = rand() % m;
    } while(indexOf(r, set, i) != -1);
    set[i] = r;

    if (i == gi_sum[j]-1){      
      weights[j] = getWeight(set, j, ks, gi_sum, kelem, n, A);      
      j++;
    }    
  }  
}


void init_soln(int *initial_soln, int *set, double *weights, int m, int n, int t, int* gi_sum, int* ks, int** kelem, int num_genes, mutation_data_t *A){
  // Declarations
  int i, j=0;

  // Select t random gene sets of size k, making sure
  // there are no duplicates
  for(i = 0; i < num_genes; i++){    
    set[i] = initial_soln[i];
    if (i == gi_sum[j]-1){      
      weights[j] = getWeight(set, j, ks, gi_sum, kelem, n, A);
      j++;
    }    
  }  
}


int get_set_to_exchange(int to_exchange, int t, int *group_index_sum){
  int i;
  if (to_exchange != -1){
    for (i=0; i < t; i++)
      if (to_exchange < group_index_sum[i])
        return i;
    return -1;
  }
  else{
    return -1;
  }

}

bool checkSubtype(int m, int next_gene, int exchanged, int set_to_exchange, int other_set_to_exchange, int *gi_sum,int *ks, int *set){
  int i, t_i, start, end;
  bool containsubtype = false;  

  // check set to exchange 
  if (next_gene >= m && exchanged < m){
    t_i = set_to_exchange;
    start = gi_sum[t_i] - ks[t_i];
    end = gi_sum[t_i];    
    for (i = start; i < end; i++) {
      if(set[i] >= m){
        containsubtype = true;
        break;
      }
    }      
  }
  if (other_set_to_exchange != -1 && other_set_to_exchange != set_to_exchange && containsubtype == false){  // swapping between gene sets
    if (next_gene < m && exchanged >= m){ // check
      t_i = other_set_to_exchange;
      start = gi_sum[t_i] - ks[t_i];
      end = gi_sum[t_i];      
      for (i = start; i < end; i++) {
        if(set[i] >= m){
          containsubtype = true;
          break;
        }
      }    
    }
  }  

  return containsubtype;
}

// Multi-Dendrix MCMC
void comet_mcmc(mutation_data_t *A, int m, int n, int *ks, int t, int num_iters, int step_len, int amp, 
                    int nt, double binom_cutoff, int *initial_soln, int size_initial_soln,
		                int size_subtype, double pval_cutoff, int **gene_sets, double **set_weights, int **set_functions, 
                    int verbose){
  // Declarations
  int i, j, step_index, next_gene, gene_index, num_genes;
  int other_set_to_exchange, set_to_exchange, to_exchange, exchanged;
  int steps, prog_step = num_iters/72; // for progress bar
  int num_swapped = 0, num_within_swap = 0; // MCMC stats
  double log_acceptance_prob, logcoin;
  double *weights, *next_weights, next_weight, current_weight = 0.0;
  int *set;
  clock_t begin = clock(), end;

  // group index for indexing irregular ks, e.g. ks = [4,3,3], group_index[[0,1,2,3],[4,5,6],[7,8,9]]  
  
  int *group_index_sum;    
  group_index_sum = malloc(sizeof(int) * t);
  for (i=0; i < t; i++){
    group_index_sum[i] = sum_array(ks, 0, i+1);     
  }
  
  // pre-compute locations of co-cells for each k
  int max_k = arr_max(ks, t);
  int **k_elem = malloc(sizeof(int *) * max_k);  

  for (i=0; i< max_k; i++){    
    k_elem[i] = get_co_cells(i+1);    
  }

  num_genes = sum_array(ks, 0, t);
  
  // Set up global variables
  mcmc_amp = amp;  
  co_cutoff = nt;
  binom_pvalthreshold = binom_cutoff;
  pvalthresh = pval_cutoff;
  permutation_iteration = (int)(1/pvalthresh);

  //printf("mcmc_amp: %d, permutation_times: %d, co-occ cutoff: %d, binom-p cutoff :%e.\n", mcmc_amp, permutation_iteration, co_cutoff, binom_pvalthreshold);
  // Allocations
  set          = malloc(sizeof(int) * num_genes);
  weights      = malloc(sizeof(weights[0]) * t);
  next_weights = malloc(sizeof(next_weights[0]) * t);

  // Intialize random number generator with the current time
  time_t T;
  srand((unsigned) time(&T));
  
  
  // Initialize the Markov chain
  if (size_initial_soln == 0){    
    rand_soln(set, weights, m, n, t, group_index_sum, ks, k_elem, num_genes, A);    
  }
  else{
    init_soln(initial_soln, set, weights, m, n, t, group_index_sum, ks, k_elem, num_genes, A);
  }

  current_weight = sum(weights, t);  
  //printf("current weight: %f, number of genes: %d\n", current_weight, num_genes);
  // for (i=0; i<num_genes;i++) printf("gene id %d: %d\n", i, set[i]);
  
  // Run MCMC process
  //c = 0.5;
  step_index = 0;  
  for (i = 0; i < num_iters; i++){
    while(1){
      next_gene        = (double) (rand() % m);
      to_exchange      = rand() % num_genes;
      set_to_exchange  = get_set_to_exchange(to_exchange, t, group_index_sum);
      gene_index       = indexOf(next_gene, set, num_genes);
      exchanged        = set[to_exchange];
      other_set_to_exchange = get_set_to_exchange(gene_index, t, group_index_sum);  
      
      if(size_subtype == 0){
        break;
      }
      // check if new sampling contains two subtype in the sample gene set 
      if (checkSubtype(m - size_subtype, next_gene, exchanged, set_to_exchange, other_set_to_exchange, group_index_sum, ks, set) == false){
          break;                  
      }        
    }
    dendrix_pass     = false;      
    // Construct the next solution
    if (gene_index != -1){ // next_gene is in current set
      set[gene_index]  = exchanged;
    }
    set[to_exchange] = next_gene;

    // Update the weights
    copyArr(weights, next_weights, t); // copy a -> b
    next_weights[set_to_exchange] = getWeight(set, set_to_exchange, ks, group_index_sum, k_elem, n, A);    
    
    if (dendrix_pass && gene_index != -1){ 
      other_set_to_exchange = get_set_to_exchange(gene_index, t, group_index_sum);
      if (set_to_exchange != other_set_to_exchange){
        next_weights[other_set_to_exchange] = getWeight(set, other_set_to_exchange, ks, group_index_sum, k_elem, n, A);
      }
    }
    if (dendrix_pass){
      
      next_weight = sum(next_weights, t);
      // Transition to the next state pseudorandomly based on the acceptance ratio
      //log_acceptance_prob = c * next_weight - c * current_weight;
      log_acceptance_prob = next_weight - current_weight;
          
      logcoin = log((double) rand() / (double) RAND_MAX);
      if (logcoin <= log_acceptance_prob){
        if (gene_index != -1) num_within_swap += 1;
        current_weight = next_weight;
        copyArr(next_weights, weights, t);
        num_swapped += 1;
      }
      else{
        dendrix_pass = false;
      }
    }
    
    if (!dendrix_pass){ 
      // dendrix < 0 or didn't pass acceptance ratio      
      if (gene_index != -1){
        set[gene_index] = next_gene;
      }
      set[to_exchange] = exchanged; // Undo the swap      
    }
  

    /* Thin out the chain by periodically saving solutions */
    if ( (i+1) % step_len == 0 ){

      for (j=0; j < num_genes; j++) gene_sets[step_index][j] = set[j];
      for (j=0; j < t; j++) {
        set_weights[step_index][j] = weights[j];        
        set_functions[step_index][j] = getFunction(set, j, ks, group_index_sum);
      }
      step_index += 1;
    }
    
    /* Output a simple progress bar */
    if ( (verbose == 1) && (i % prog_step == 0) ){
      steps = (int) 72. * i / num_iters;
      for (j = 0; j < steps; j++) printf("+");
      for (j = 0; j < 72 - steps; j++) printf(" ");
      printf(" [ %d%% ]\r", (int) (100. * i / num_iters));
      fflush(stdout);
    }
  }
  if (verbose == 1){
    printf("\n");
    printf("  - Swapped percentage: %d / %d (%d within)\n", num_swapped,
           num_iters, num_within_swap);
    printf("  - Look ups: %d\n", lookups);
    printf("  - Unique gene sets: %d\n", calculations);
    end = clock();
    printf("  - MCMC runtime: %f secs\n", (double)(end - begin) / CLOCKS_PER_SEC);
  }

  // Free memory
  free(set);
  free(weights);
  free(next_weights);
  free(group_index_sum);
  for (i = 0; i < max_k; i++)
    if(k_elem[i] != NULL) free(k_elem[i]);
  
  
  //delete_all();  // free uthash

}


