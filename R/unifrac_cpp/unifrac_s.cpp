/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "tree_s.hpp"
#include "biom_interface_s.hpp"
#include "unifrac_s.hpp"
#include <unordered_map>
#include <cstdlib>
#include <thread>
#include <signal.h>
#include <stdarg.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>

#include "unifrac_internal_s.hpp"

#include <Rcpp.h>

using namespace su;

// Computes Faith's PD for the samples in  `table` over the phylogenetic
// tree given by `tree`.
// Assure that tree does not contain ids that are not in table
std::vector<double> su::faith_pd(tse_interface &table,
                  BPTree &tree) {
    PropStack<double> propstack(table.n_samples); // construction seems to go okay

    
    uint32_t node;
    std::vector<double> node_proportions;
    double length;
    
    std::vector<double> results = std::vector<double>();

    // for node in postorderselect
    for(unsigned int k = 0; k < (tree.nparens / 2) - 1; k++) {
        node = tree.postorderselect(k);
        
        // get branch length
        length = tree.lengths[node];
        
        // get node proportions and set intermediate scores
        
        node_proportions = set_proportions(tree, node, table, propstack);
        
        for (unsigned int sample = 0; sample < table.n_samples; sample++){
            // calculate contribution of node to score
            results.push_back((node_proportions[sample] > 0) * length);
        }
    }
    return results;
    //return(std::vector<double>());
}
