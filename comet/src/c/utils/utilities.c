#include "../cometmodule.h"

/*******************************************************************************
  Utilities 
 *******************************************************************************/

double lnfac(int a){
    double z;
    int y;
    if ((a == 1) || a == 0) return 0.0;
    else{
        z = 0;
        for (y = 2; y<=a; y++ ) z = log(y)+z;
        return z;
    }
}

// Sort a pair of integers in ascending order
int ascending(const void * elem1, const void * elem2){
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

void precompute_factorials(int N){
    int i;
    lnfacs = malloc(sizeof(double) * (N+1));
    for (i=0; i < N+1; i++) {
		lnfacs[i] = lnfac(i);  
	}
}

/* Returns the index of n in array where k is the length of array, or -1 if n is
   not in array.*/
int indexOf(int n, int *array, int k){
    int i;
    for (i=0; i < k; i++) if(array[i] == n) return i;
    return -1;
}

// Numeric functions
int min(int x, int y){ return x < y ? x : y; }
int max(int x, int y){ return x > y ? x : y; }
int arr_min(int *arr, int len) {
    int i, min;
    if (arr == NULL || len < 1) {
        return 0;
    }
    min = arr[0];
    for (i=0; i<len; i++) {
        min = arr[i] < min ? arr[i] : min;        
    } 
    return min;
}
int arr_max(int *arr, int len) {
    int i, max;
    if (arr == NULL || len < 1) {
        return 0;
    }
    max = arr[0];
    for (i=0; i<len; i++) {
        max = arr[i] > max ? arr[i] : max;        
    } 
    return max;
}

/* Comparison function for sorting from greatest to smallest. */
/*int comp_g2s_numMutations(const void *e1, const void *e2){
    int f = glbl_G2numMuts[ *((int *)e1) ];
    int s = glbl_G2numMuts[ *((int *)e2) ];
    if (f > s) return -1; // sort greatest to smallest
    if (f < s) return 1;
    return 0;
}*/

int sum_int_array(int *arr, int len) {
    int i, sum;
    sum = 0;
    for (i=0; i<len; i++) {
        sum += arr[i];
    }
    return sum;
}

// Array functions
/* Converts a 2D array into a PyObject of type PyList */
void convert_2Darr_pylist(int rows, int cols, PyObject *pyArr, int **arr) {
    int i;
    PyObject *this_row;
    for (i=0; i < rows; i++){
        this_row = PyList_New(cols);
        convert_int_array_to_pylist(cols, this_row, arr[i]);
        PyList_SET_ITEM(pyArr, i, this_row);
    }
}

/* Converts an array into a PyObject of type PyList */
void convert_int_array_to_pylist(int len, PyObject *pyArr, int *arr) {
    int j;
    if (arr == NULL) {
        arr = malloc(sizeof(int) * len);
        memset(arr, 0, sizeof(int) * len);
    }
    for (j=0; j < len; j++) {
        PyList_SET_ITEM(pyArr, j, PyInt_FromLong(arr[j]));
    }
}
/* Converts an array into a PyObject of type PyList */
void convert_double_array_to_pylist(int len, PyObject *pyArr, double *arr) {
    int i;
    if (arr == NULL) {
        arr = malloc(sizeof(double) * len);
        memset(arr, 0, sizeof(double) * len);
    }

    for (i=0; i<len; i++) {
        PyList_SET_ITEM(pyArr, i, PyFloat_FromDouble(arr[i]));
    }
}

/* Save array arr at index i of frozenarrays. Array is length len.*/
void freeze_int_array(int i, int len, int **frozenarrays, int *arr) {
    frozenarrays[i]  = malloc(sizeof(int) * len);
    memcpy(frozenarrays[i], arr, (sizeof(int) * len));
}

/* Free the pointers in the array, then free the array. */
void free_ptr_array(void **array, int len) {
    int i;
    for (i=0; i < len; i++) {
        free(array[i]);
    }
    free(array);
}

/* Convert all frozen sets to PyObject to return to python program*/
void frozensets2pylists(frozen_arrays_t *frozen, int numFrozen, int k, PyObject *solns, PyObject *weights, PyObject *tables, PyObject *probs) {
    /* Store the frozensets in a Python list for returning */
    convert_2Darr_pylist(numFrozen, k, solns, frozen->frozensets);
    convert_double_array_to_pylist(numFrozen, weights, frozen->frozenadjweights);
    convert_2Darr_pylist(numFrozen, 1 << k, tables, frozen->frozentables);
    convert_double_array_to_pylist(numFrozen, probs, frozen->frozenprobs);
}


/* Allocate all memory needed for frozen value arrays for dendrix */
frozen_arrays_t *frozen_allocate(int numFrozen) {
    frozen_arrays_t *frozen = (frozen_arrays_t*) malloc(sizeof(frozen_arrays_t));
    frozen->frozensets    = malloc(numFrozen * sizeof(int *));
    frozen->frozenadjweights = malloc(numFrozen * sizeof(double));
    frozen->frozentables  = malloc(numFrozen * sizeof(int *));
    frozen->frozenprobs   = malloc(numFrozen * sizeof(double));
    return frozen;
}

/* Free all memory associated with frozen value arrays for dendrix */
void frozen_free(int numFrozen, frozen_arrays_t *frozen) {
    free_ptr_array((void **) frozen->frozensets, numFrozen);
    free(frozen->frozenadjweights);
    free_ptr_array((void **) frozen->frozentables, numFrozen);
    free(frozen->frozenprobs);
    free(frozen);
}

current_soln_values_t *current_arrays_allocate(int k) {
    current_soln_values_t *curr_soln_vals = (current_soln_values_t *) malloc(sizeof(current_soln_values_t));
    curr_soln_vals->current_soln    = malloc(sizeof(int) * k);
    curr_soln_vals->current_table   = (int *) malloc(sizeof(int) * 1 << k);
    curr_soln_vals->current_prob = 0.0;
    curr_soln_vals->current_weight = 0.0;
    curr_soln_vals->current_num_tbls = 0;
    return curr_soln_vals;
}

void current_arrays_free(current_soln_values_t *curr_soln_vals) {
    free(curr_soln_vals->current_soln);
    free(curr_soln_vals->current_table);
    free(curr_soln_vals);
}

/*******************************************************************************
 * Memoization
 *******************************************************************************/
// Functions for specifying the behavior of solution-probability hash table 
/* For a set of genes, creates a string representation of the set stored in hash */
/*
void set_to_string(int *array, int len, char *hash) {
   int *tmp_array;
   int a, b;
   int comp(const void *e1, const void *e2){
       a = *((int *)e1);
       b = *((int *)e2);
       if (a > b) return -1; // sort greatest to smallest
       if (a < b) return 1;
       return 0;
   }
   tmp_array = (int *) malloc(sizeof(int) * len);
   memcpy(tmp_array, array, sizeof(int) * len);
   qsort(tmp_array, len, sizeof(*tmp_array), comp);
   int_arr_to_string(tmp_array, len, hash, " "); 
   free(tmp_array);
}
*/
void int_arr_to_string(int *array, int len, char *dest, char *delimiter) {
    int i;
    int dest_len = len + len - 1; //len numbers separated by spaces
    char *tmp = (char *) malloc(sizeof(char) * 10); //Assume each gene number at most 10 digits
    
    for (i=0; i<dest_len; i++) {
        if (i%2 == 0) {
            snprintf(tmp, 100, "%d", array[(i+1)/2]);
            strncat(dest, tmp, 100);
        } else {
            strncat(dest, delimiter, 100); 
        }
    }
    free(tmp);
}

// Copy an array of doubles
void copyArr(double *A, double *B, int n){
  int i;
  for (i = 0; i < n; i++) B[i] = A[i];
}

// Sum an array of doubles
double sum(double *arr, int n){
  int i;
  double total = 0;
  for (i = 0; i < n; i++) total += arr[i];
  return total;
}

int sum_array(int *arr, int start, int end){
  int i, sum=0;
  for (i=start; i<end; i++){
    sum += arr[i];
    //printf("%d\n", arr[i]);
  }
  return sum;
}

