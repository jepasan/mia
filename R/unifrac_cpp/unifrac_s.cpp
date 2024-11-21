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

using namespace su;

std::string su::test_table_ids_are_subset_of_tree(su::biom_interface &table, su::BPTree &tree) {
    std::unordered_set<std::string> tip_names = tree.get_tip_names();
    std::unordered_set<std::string>::const_iterator hit;
    std::string a_missing_name = "";

    for(auto i : table.obs_ids) {
        hit = tip_names.find(i);
        if(hit == tip_names.end()) {
            a_missing_name = i;
            break;
        }
    }

    return a_missing_name;
}



// Computes Faith's PD for the samples in  `table` over the phylogenetic
// tree given by `tree`.
// Assure that tree does not contain ids that are not in table
void su::faith_pd(biom_interface &table,
                  BPTree &tree,
                  double* result) {
    PropStack<double> propstack(table.n_samples);

    uint32_t node;
    double *node_proportions;
    double length;

    // for node in postorderselect
    for(unsigned int k = 0; k < (tree.nparens / 2) - 1; k++) {
        node = tree.postorderselect(k);
        // get branch length
        length = tree.lengths[node];

        // get node proportions and set intermediate scores
        node_proportions = propstack.pop(node);
        set_proportions(node_proportions, tree, node, table, propstack);

        for (unsigned int sample = 0; sample < table.n_samples; sample++){
            // calculate contribution of node to score
            result[sample] += (node_proportions[sample] > 0) * length;
        }
    }
}
