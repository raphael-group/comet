#include "cometmodule.h"

char weight_func_c;

double calc_soln_weight(mutation_data_t *mut_data, current_soln_values_t *curr_arrays, int k, double pvalthresh);
void enumerate_helper(mutation_data_t *mut_data, int k, current_soln_values_t *curr_soln_vals,
                      frozen_arrays_t *frozen, int start_i, int *step,
                      int layer, double pvalthresh); 
int enumerate(mutation_data_t *mut_data, int k, frozen_arrays_t *frozen, double pvalthresh);

/* Calculate the weight of a gene set */
double calc_soln_weight(mutation_data_t *mut_data, current_soln_values_t *curr_arrays,
                        int k, double pvalthresh) {
    int ctbl_len = 1 << k; 
    int final_num_tbls;
    int *tmp_soln = malloc(sizeof(int) * k);
    int *ctbl = malloc(sizeof(int) * ctbl_len); 

    //glbl_G2numMuts = mut_data->G2numMuts; /* Set global G2numMuts array for sorting purposes 

    memcpy(tmp_soln, curr_arrays->current_soln, sizeof(int) * k);

    /* Sort the solution by the number of mutations */
    //qsort(tmp_soln, k, sizeof(*tmp_soln), comp_g2s_numMutations);
    contingency_tbl_gen(mut_data, k, ctbl, tmp_soln);

    memcpy(curr_arrays->current_table, ctbl, (sizeof(int) * ctbl_len));
    //curr_arrays->current_weight = weight_func(mut_data, curr_arrays, ctbl, k, pvalthresh);
    if (weight_func_c == 'e'){
        curr_arrays->current_prob = comet_exact_test(k, mut_data->numPatients, ctbl, &final_num_tbls, pvalthresh);
    }
    else if (weight_func_c == 'd'){
        curr_arrays->current_prob = dendrix(ctbl, k);
    }
    else if (weight_func_c == 'b'){
        curr_arrays->current_prob = comet_binomial_test(k, mut_data->numPatients, ctbl, pvalthresh);   
    }
    else if (weight_func_c == 'c'){
        curr_arrays->current_prob = comet_binomial_co_test(k, mut_data->numPatients, ctbl, pvalthresh);      
    }

    free(ctbl);
    free(tmp_soln);
    return 0-log(curr_arrays->current_prob);
} 

/*******************************************************************************
* Exhaustive enumeration of all gene sets of size k
*******************************************************************************/
int choose(int N, int k) {
    return ceil(exp(lnfacs[N] - (lnfacs[k] + lnfacs[N-k])));
}

int enumerate(mutation_data_t *mut_data, int k, frozen_arrays_t *frozen, double pvalthresh) {
    int step = 0;
    current_soln_values_t *curr_soln_vals;
    curr_soln_vals = current_arrays_allocate(k);
    enumerate_helper(mut_data, k, curr_soln_vals, frozen, 0, &step, 0, pvalthresh);
    current_arrays_free(curr_soln_vals); 
    return step;
}

void enumerate_helper(mutation_data_t *mut_data, int k, current_soln_values_t *curr_soln_vals,
                      frozen_arrays_t *frozen, int start_i, int *step, int layer,
                      double pvalthresh) {
    int i;
    if (layer == k) {
        double weight;
        weight = calc_soln_weight(mut_data, curr_soln_vals, k, pvalthresh);
        //printf("weight: %f", weight);
        freeze_int_array(*step, 1 << k, frozen->frozentables, curr_soln_vals->current_table);
        freeze_int_array(*step, k, frozen->frozensets, curr_soln_vals->current_soln);
        frozen->frozenprobs[*step] = curr_soln_vals->current_prob;
        frozen->frozenadjweights[*step] = weight;
        *step = *step + 1;
        return;
    } else {
        for (i=start_i; i < mut_data->numGenes; i++) {
            curr_soln_vals->current_soln[layer] = i;
            enumerate_helper(mut_data, k, curr_soln_vals, frozen, i+1, step, layer+1,
                             pvalthresh);
        }
    }
}

frozen_arrays_t *comet_exhaustive(mutation_data_t *mut_data, int k, int *numSets, double pvalthresh) {
    //const char *logfile_name = "./exhaustive_log.tsv";
    frozen_arrays_t *frozen;
    frozen = frozen_allocate(*numSets);
    //if (log_results > 0) {
    //    init_io_globals(logfile_name);
    //}
    /* Fix potential off-by-one from the "choose" function */    
    *numSets = enumerate(mut_data, k, frozen, pvalthresh); 
    //if (log_results > 0) {
    //    fclose(logfile);
    //}
    return frozen;
}
