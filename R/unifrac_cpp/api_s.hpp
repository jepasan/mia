#include "task_parameters.hpp"

#ifdef __cplusplus
#include <vector>
#include <Rcpp.h>
#define EXTERN extern "C"

#else
#include <stdbool.h>
#define EXTERN
#endif

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
std::vector<double> faith_pd_one_off(const Rcpp::S4 & treeSE, bool rooted);