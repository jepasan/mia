/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include <stack>
#include <vector>
#include <unordered_map>
#include <thread>
#include <pthread.h>

#ifndef __UNIFRAC

#include "task_parameters.hpp"
#include "biom_interface.hpp"

    namespace su {
    
        void faith_pd(biom_interface &table, BPTree &tree, double* result);
    
        std::string test_table_ids_are_subset_of_tree(biom_interface &table, BPTree &tree);
        
    }
    
#define __UNIFRAC 1
#endif
