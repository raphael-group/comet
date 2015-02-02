#include "cometmodule.h"
#include "utils/cephes/mconf.h"
#include <time.h>

// Contingency table functions shared by exact and binomial tests
//int sum_cells( int *tbl, int *cells, int num_entries );
void fixed_cells( int k, int pos, int val, int *cells );
double denom(int k, int num_entries, int N, int *tbl);
int* random_sample(int max, int x);

// Exact test functions
void derive_remaining_cells(int k, int N, int *margins, int *ex_cells, int *tbl, int *mar_rems);
int min_affected_margin(int k, int cell, int *mar_rems);
int exact_test_helper(double *pval, int *num_tbls, int k, double pvalthresh, int num_entries,
                      int N, double numerator, int *margins, int *ex_cells, int *co_cells,
                      int num_co_cells, int *tbl, int **mar_stack, int co_in, int T_rem, int T_obs);

// Constants for 
static double const OVER_THRESH  = -1;
static double const UNDER_THRESH = 0;
//static double const OVER_WEIGHT  = 0.04;


/*******************************************************************************
  Permutation functions
 *******************************************************************************/

int* random_sample(int max, int x)
{   
   int n = x;
   int N = max+1;
   int *a;
   a = malloc(x*sizeof(int));

   int t = 0; // total input records dealt with
   int m = 0; // number of items selected so far
   double u;

   while (m < n)
   {
       u = rand() / (RAND_MAX + 1.0) ;
       if ( (N - t)*u >= n - m )
       {
          t++;
       }
       else
       {
          a[m] = t;
          t++; m++;
       }
   }
   return a;
}

/*******************************************************************************
  Weight functions
 *******************************************************************************/


double comet_permutation_test(int k, int num_permutation, int num_patients, int *gene_set, int *freq, int *tbl){
  // Parameters
  int i, p2g_index, l; // Helper variables  

  int *ex_cells     = get_ex_cells(k);
  int Tobs = sum_cells(tbl, ex_cells, k); 
  free(ex_cells);      
   
  int count = 0;
  int count_mid_p = 0;
    
  for (i=0; i<num_permutation; i++){
    // initialization
    int Trand = 0;
    int permute_mut_data[num_patients];    
    for(l=0; l<num_patients; l++){
      permute_mut_data[l] = 0;
    }
    // randomize mutation data in permute_mut_data  
    for(p2g_index=0; p2g_index<k; p2g_index++){        
        // random sampling from a list with freq samples?                
        int *new_samples = random_sample(num_patients-1, freq[p2g_index]);        
        for (l=0; l<freq[p2g_index]; l++){                    
          permute_mut_data[new_samples[l]]++;          
        }                
        free(new_samples);
    }
    for(l=0; l<num_patients; l++){
      if (permute_mut_data[l] == 1) Trand++;
    }    

    if (Tobs <= Trand){
      count++;
    }
    if (Tobs < Trand){
      count_mid_p++;
    }
  }  
  
  return (double)(count+count_mid_p)/2;

}


/* General exact test for Dendrix ++*/
int sum_cells( int *tbl, int *cells, int num_entries ){
  int i, sum = 0;
  for (i = 0; i < num_entries; i++){
    sum += tbl[cells[i]];
  }
  return sum;
}

void fixed_cells( int k, int pos, int val, int *cells ){
  // Increment pos so counting starts at one
  pos += 1;
  
  // Define variables
  int *before, *after;
  int i, j, cell_num;

  int num_before  = pow(2, pos-1);
  int num_after   = pow(2, k-pos);
  
  // Initialize arrays
  before = malloc(sizeof(int) * num_before);
  after = malloc(sizeof(int) * num_after);
  
  // Compute all integers of k bits before and after the fixed position
  for (i=0; i < num_before; i++) before[i] = i;
  for (i=0; i < num_after; i++) after[i]  = i << pos;	
  
  // Combine the binary strings and set the bit in the given position to 1
  cell_num = 0;
  for (i = 0; i < num_before; i++){
    for (j=0; j < num_after; j++){
      cells[cell_num] = (before[i] + after[j]) | (val << (pos-1));
      cell_num += 1;
    }
  }
  
  // Clean up memory
  free(before);
  free(after);

}


double denom(int k, int num_entries, int N, int *tbl){
  int i;
  double total = (k-1)*lnfacs[N];
  for (i=0; i < num_entries; i++){ total += lnfacs[tbl[i]]; }
  return total;
}

// Create a list of all k-bit binary strings with one 1
int *get_ex_cells(int k){
  int *ex_cells;
  int i;
  ex_cells = malloc(sizeof(int) * k);
  for (i=0; i < k; i++) ex_cells[i] = 1 << i;
  return ex_cells;
}

// Counts the number of ones in the binary representation of an integer
int num_ones(int n){
  int count = 0;
  while (n > 0){
    count += 1;
    n &= n - 1;
  }
  return count;
}

// Create a list of all k-bit binary strings with >1 1s
int *get_co_cells(int k){
  int *co_cells;
  int i, cell_num, num_co_cells = pow(2, k)-k-1;
  co_cells = malloc(sizeof(int) * num_co_cells);
  cell_num = 0;
  for (i=pow(2, k); i > 0; i--){
    if( num_ones(i) > 1){
      co_cells[cell_num] = i;
      cell_num++;
    }
  }
  return co_cells;
}

void derive_remaining_cells(int k, int N, int *margins, int *ex_cells, int *tbl, int *mar_rems){
  int i;
  for (i = 0; i < k; i++){
    tbl[ex_cells[i]] = mar_rems[i];
  }
  tbl[0] = N;
  for (i = 1; i < pow(2, k); i++) tbl[0] -= tbl[i];
}

int min_affected_margin(int k, int cell, int *mar_rems){
  int i, min_val;
  min_val = INT_MAX;
  for (i = 0; i < k; i++){
    // Check if the ith variable is one for the given cell
    // and replace if it is less than the current min
    if ( cell & (1 << i) && mar_rems[i] <= min_val ) min_val = mar_rems[i];
  }
  return min_val;
}

// Determines if any integer in an array is negative
int contains_negative(int *arr, int len){
  int i;
  for (i = 0; i < len; i++){
    if (arr[i] < 0) return 1;
  }
  return 0;
}

int exact_test_helper(double *pval, int *num_tbls, int k, double pvalthresh, int num_entries,
                      int N, double numerator, int *margins, int *ex_cells, int *co_cells,
                      int num_co_cells, int *tbl, int **mar_stack, int co_in, int T_rem, int T_obs){
  int res = UNDER_THRESH;
  if (co_in >= num_co_cells){
    derive_remaining_cells( k, N, margins, ex_cells, tbl, mar_stack[co_in] );
    if (contains_negative(tbl, num_entries) == 0){
      double addp = exp( numerator - denom( k, num_entries, N,  tbl) );
      pval[0] += addp;
      if (T_obs < sum_cells(tbl, ex_cells, k)){ // T > T_x
        pval[1] += addp;
      }
      num_tbls[0] += 1;
    }
  if ((pval[0]+pval[1])/2 > pvalthresh) {	
      res = OVER_THRESH;
    }
  }
  else {
    // Define required variables
    int i, cell, val, MarRem;
    double coef;
    int *mar_rems;
    
    cell     = co_cells[co_in];
    coef     = num_ones( cell );
    mar_rems = mar_stack[co_in];
    
    // Determine which variables are in the margin
    MarRem = min_affected_margin( k, cell, mar_rems );
    
    // Iterate over the possible values the current cell can take
    for (val = 0; val < min(MarRem, (int) floor(T_rem/coef)) + 1; val++){
      // Update margins
      for (i=0; i < k; i++){
        if (cell & (1 << i)) mar_stack[co_in+1][i] = mar_rems[i] - val;
        else mar_stack[co_in+1][i] = mar_rems[i];
      }
      
      // Create new table using the current value
      tbl[cell] = val;
      res = exact_test_helper( pval, num_tbls, k, pvalthresh, num_entries, N, numerator,
                               margins, ex_cells, co_cells, num_co_cells, tbl,
                               mar_stack, co_in + 1, T_rem-coef*val, T_obs);
      if (res < 0) {
        break;
      }
    }
  }
  return res;
}

double comet_exact_test(int k, int N, int *ctbl, int *final_num_tbls, double pvalthresh){
  // Computation variables
  int kbar, Tobs, T, margin_size, num_co_cells;
  double numerator;
  double *pval;
  int *ex_cells, *co_cells, *blank_tbl, *margins, *num_tbls, *cells;
  int **mar_stack;
  
  int num_entries, i, res; // Helper variables
  double final_p_value;
  
  num_entries = 1 << k;
  
  // Compute binary representation of each cell
  ex_cells     = get_ex_cells(k);
  co_cells     = get_co_cells(k);
  num_co_cells = pow(2, k)-k-1;
  
  // Allocate memory stack for marginal remainders
  mar_stack = malloc(sizeof(int *) * (num_co_cells + 1));
  for (i=0; i < num_co_cells + 1; i++) mar_stack[i] = malloc(sizeof(int) * k);
  cells = malloc(sizeof(int) * 1 << (k-1));
  
  /* Compute additional fixed values */
  // Margins
  margins = malloc(sizeof(int) * 2 * k);
  margin_size = pow(2, k-1);
  for (i = 0; i < k; i++){
    // Negative margin
    fixed_cells(k, i, 0, cells);
    margins[i]   = sum_cells( ctbl, cells, margin_size );
    
    // Positive margin
    fixed_cells(k, i, 1, cells);
    margins[i+k] = sum_cells( ctbl, cells, margin_size );
    
    // Initialize margin stack
    mar_stack[0][i] = margins[i+k];
  }
  
  // Numerator
  numerator = 0.0;
  for (i = 0; i < 2*k; i++) numerator += lnfacs[margins[i]];
  
  // Observed
  Tobs = sum_cells(ctbl, ex_cells, k);
  kbar = 0;
  for (i = 0; i < k; i++) kbar += margins[i+k];
  
  /* Set up recursion */
  // Set remaining co-occurrences allowed
  T = kbar - Tobs;
  pval = malloc(sizeof(double));
  num_tbls = malloc(sizeof(int));
  pval[0] = 0.0;
  pval[1] = 0.0; // mid-pvalue
  num_tbls[0] = 0;
  
  // Construct a blank table to start with
  blank_tbl = malloc(sizeof(int) * num_entries);
  
  // Run the exact test recursion which will update the results array
  // with the number of extreme tables and the pval
  res = exact_test_helper( pval, num_tbls, k, pvalthresh, num_entries, N, numerator,
                           margins, ex_cells, co_cells, num_co_cells, blank_tbl,
                           mar_stack, 0, T, Tobs);
  
  *final_num_tbls = num_tbls[0];
  final_p_value = (pval[0]+pval[1])/2;
  
  // Free memory
  free(cells);
  free(margins);
  free(pval);
  free(num_tbls);
  free(co_cells);
  free(ex_cells);
  free(blank_tbl);
  free_ptr_array((void **) mar_stack, num_co_cells + 1);
  
  return (res == OVER_THRESH) ? res : final_p_value;

}

/*******************************************************************************
  Dendrix++ binomial approximation
 *******************************************************************************/


double binomial_cdf(k, n, p)
int k, n;
double p;
{
    double dk, dn;
    
    if (k == n) return (1.0);
    dn = n - k;
    if (k == 0) dk = pow(1.0 - p, dn);
    else {
      dk = k + 1;
      dk = incbet(dn, dk, 1.0 - p);
    }
    return (dk);
}


// C version of Dendrix++ binomial approximation
double comet_binomial_test(int k, int N, int *tbl, double pvalthresh){
  // Declarationa
  int i, j, margin_size, Tobs;
  int *cells, *ex_cells;
  double p_e, prod, pval, num_patients, pval_m;
  double *margin_probs;

  /* Compute additional fixed values */
  // Margins
  margin_probs = malloc(sizeof(double) * 2 * k);
  margin_size = 1 << (k-1);
  num_patients = (double) N;
  cells = malloc(sizeof(int) * 1 << (k-1));
  for (i = 0; i < k; i++){
    // Negative margin
    fixed_cells(k, i, 0, cells);
    margin_probs[i] = sum_cells( tbl, cells, margin_size ) / num_patients;
    
    // Positive margin
    fixed_cells(k, i, 1, cells);
    margin_probs[i+k] = sum_cells( tbl, cells, margin_size ) / num_patients;
  }
  
  // Exclusive cells
  ex_cells = get_ex_cells(k);
  
  // Observed
  Tobs = sum_cells(tbl, ex_cells, k);
  
  // Compute the probability of an exclusive mutation: p_e
  p_e = 0.0;
  for (i=0; i < k; i++){
    prod = 1.0;
    for (j=0; j < k; j++){
      if (i==j) prod *= margin_probs[j+k];
      else  prod *= 1-margin_probs[j+k];
    }
    p_e += prod;
  }
  
  // Tobs ~ Binomial(N, p_e) => Pr(T > Tobs) = 1-BinomCDF(Tobs; N, p_e)
  //pval = 1.0 - binomial_cdf( Tobs-1, N, p_e );
  pval = 1.0 - binomial_cdf(Tobs-1, N, p_e);
  pval_m = 1.0 - binomial_cdf(Tobs, N, p_e);
  
  // Free memory
  free(cells);
  free(margin_probs);
  free(ex_cells);
  
  return 0.5*(pval+pval_m);
}

double comet_binomial_co_test(int k, int N, int *tbl, double pvalthresh){
  // Declarationa
  int i, margin_size, Obs;
  int *cells;
  double p_c, pval, num_patients;
  double *margin_probs;

  /* Compute additional fixed values */
  // Margins
  margin_probs = malloc(sizeof(double) * 2 * k);
  margin_size = 1 << (k-1);
  num_patients = (double) N;
  cells = malloc(sizeof(int) * 1 << (k-1));
  for (i = 0; i < k; i++){
    // Negative margin
    fixed_cells(k, i, 0, cells);
    margin_probs[i] = sum_cells( tbl, cells, margin_size ) / num_patients;
    
    // Positive margin
    fixed_cells(k, i, 1, cells);
    margin_probs[i+k] = sum_cells( tbl, cells, margin_size ) / num_patients;
  }
  
  // Observed: number of co-occurring
  Obs = tbl[(1 << k) -1];
  
  // Compute the probability of a completely co-occurring mutation
  p_c = 1.0;
  for (i = 0; i < k; i++) p_c *= margin_probs[i+k];
  
  // Tobs ~ Binomial(N, p_e) => Pr(T > Tobs) = 1-BinomCDF(Tobs; N, p_e)
  pval = 1.0 - binomial_cdf( Obs-1, N, p_c );
  
  // Free memory
  free(cells);
  free(margin_probs);
  
  return pval;
}


void comet_phi(int *genes, int k, int n, mutation_data_t *A, int *ctbl, int** kelem, double *score, int *func, int co_cutoff, int permutation_iteration, double binom_pvalthresh, double exact_pvalthresh){
  /* 
    key->genes: genes in the set, for permutation test
    k: the number of genes
    n: the number of total samples
    A: mutation data, for permutation test
    ctbl: contingency table, for exact and binom test
    kelem: index of co-occurring cells in 2^k ctbl
    &score: output phi
    &func: index for weight function that is used for calculating the phi 
    (1: exact, 2: binom, 3: permutation, 4: exact with fewer 100 tables)
    co_cutoff: cutoff for co-occurring (determine binom test will be a good approximation for exact test)
    permutation_iteration: number of permutation test
    binom_pvalthresh: cutoff for binom test
    exact_pvalthresh: cutoff for exact test (stop table enumeration if sum of table probabilities is larger than exac_pvalthresh)
  */  
  double permute_count = 0.0;
  double binom_pval, exact_pval;
  int final_num_tbls; // store number of enumerating tables
  int i, conum=0; // num of cooccurring mutations
  int freq[k]; // for permutation test
  int numTableCutoff = 0;

  if (k > 3){ // heuristic pipeline when k > 3    
    
    binom_pval = comet_binomial_test(k, n, ctbl, 1.1);    
    conum = sum_cells(ctbl, kelem[k-1], pow(2, k)-k-1 );   // num of cooccurring mutations 
    if (conum  > co_cutoff || binom_pval > binom_pvalthresh){         
        *score =binom_pval;
        *func = 2;
    }
    else{     
      exact_pval = comet_exact_test(k, n, ctbl, &final_num_tbls, exact_pvalthresh);
      if (exact_pval == -1){ // stop when pval > 0.01, do permutation test
    
        for (i=0; i<k; i++){
          freq[i] = A->G2numMuts[genes[i]]; // assign freq for k genes                   
        }
              
        permute_count = comet_permutation_test(k, permutation_iteration, n, genes, freq, ctbl);        
        if  (permute_count == 0.0 ){
          *score = binom_pval;
          *func = 2;
        } 
        else{ 
          *score = permute_count/permutation_iteration;
          *func = 3;
        }
      }  
      else{
        *score = exact_pval;
        if (final_num_tbls > numTableCutoff){
          *func = 1;
        }
        else{
          *func = 4;
        }
      }
    }
  } else{    
    *score = comet_exact_test(k, n, ctbl, &final_num_tbls, 1.1);
    if (final_num_tbls > numTableCutoff){
      *func = 1;
    }
    else{
      *func = 4;
    }
  }
}


/* The dendrix weight function. */
double dendrix(int *ctbl, int k) {
  int i, num_o;
  int coverage_set, sum_coverage_genes;
  int weight;
  int len = 1 << k;
  coverage_set = sum_int_array(ctbl, len) - ctbl[0];
  sum_coverage_genes = 0; 
  for (i=1; i < len; i++) {
    num_o = num_ones(i);
    sum_coverage_genes += num_o * ctbl[i];
  }
  
  weight = 2*coverage_set - sum_coverage_genes;
  return weight;
}
