#include "api_s.hpp"
#include "biom_s.hpp"
#include "tree_s.hpp"
#include "unifrac_s.hpp"
#include <fstream>
#include <iomanip>
#include <thread>
#include <cstring>
#include <stdlib.h> 

#include <fcntl.h>
#include <unistd.h>

#include <Rcpp.h>

using namespace su;
using namespace std;

// https://stackoverflow.com/a/19841704/19741
bool is_file_exists(const char *fileName) {
    std::ifstream infile(fileName);
        return infile.good();
}

void initialize_results_vec(r_vec* &result, biom& table){
    // Stores results for Faith PD
    result = (r_vec*)malloc(sizeof(results_vec));
    result->n_samples = table.n_samples;
    result->values = (double*)malloc(sizeof(double) * result->n_samples);
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);

    for(unsigned int i = 0; i < result->n_samples; i++) {
        size_t len = table.sample_ids[i].length();
        result->sample_ids[i] = (char*)malloc(sizeof(char) * len + 1);
        table.sample_ids[i].copy(result->sample_ids[i], len);
        result->sample_ids[i][len] = '\0';
        result->values[i] = 0;
    }

}

void destroy_results_vec(r_vec** result) {
    // for Faith PD
    for(unsigned int i = 0; i < (*result)->n_samples; i++) {
        free((*result)->sample_ids[i]);
    };
    free((*result)->sample_ids);
    free((*result)->values);
    free(*result);
}


/*
#define PARSE_SYNC_TREE_TABLE(tree_filename, table_filename) std::ifstream ifs(tree_filename);                                        \
 std::string content = std::string(std::istreambuf_iterator<char>(ifs),                                                               \
 std::istreambuf_iterator<char>());                                                                                                   \
 su::BPTree tree = su::BPTree(content);                                                                                               \
 su::biom table = su::biom(biom_filename);                                                                                            \
 if(table.n_samples <= 0 | table.n_obs <= 0) {                                                                                        \
 return table_empty;                                                                                                                  \
 }                                                                                                                                    \
 std::string bad_id = su::test_table_ids_are_subset_of_tree(table, tree);                                                             \
 if(bad_id != "") {                                                                                                                   \
 return table_and_tree_do_not_overlap;                                                                                                \
 }                                                                                                                                    \
 std::unordered_set<std::string> to_keep(table.obs_ids.begin(),                                                                       \
 table.obs_ids.end());                                                                                                                \
 su::BPTree tree_sheared = tree.shear(to_keep).collapse();                                                                            \
 */

compute_status faith_pd_one_off(const Rcpp::S4 & treeSE, r_vec** result, std::string newick){
    
    // Check that tree and table are non-empty and match before calling the c++ code
    // shear the tree (to contain only the obs in the table?) - Also should be done before the call?
    
    su::BPTree tree = su::BPTree(newick);
    
    initialize_results_vec(*result, table);

    // compute faithpd
    su::faith_pd(table, tree_sheared, std::ref((*result)->values));

    return okay;
}