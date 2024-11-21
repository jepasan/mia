#include "task_parameters.hpp"

#ifdef __cplusplus
#include <vector>
#include <Rcpp.h>
#define EXTERN extern "C"

#else
#include <stdbool.h>
#define EXTERN
#endif

typedef enum compute_status {okay=0, tree_missing, table_missing, table_empty, unknown_method, table_and_tree_do_not_overlap, output_error} ComputeStatus;

/* a result vector
 *
 * n_samples <uint> the number of samples.
 * values <double*> the score values of length n_samples.
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct results_vec{
    unsigned int n_samples;
    double* values;
    char** sample_ids;
} r_vec;

void destroy_results_vec(r_vec** result);

/* compute Faith PD
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * result <r_vec**> the resulting vector of computed Faith PD values
 *
 * faith_pd_one_off returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus faith_pd_one_off(const Rcpp::S4 & treeSE,
                                      r_vec** result, std::string newick);
