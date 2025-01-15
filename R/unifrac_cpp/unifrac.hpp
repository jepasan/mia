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

#include <Rcpp.h>

#ifndef __UNIFRAC

#include "biom_interface.hpp"
#include "tree.hpp"

    namespace su {
        std::vector<double> faith_pd(tse_interface &table, su::BPTree &tree);
    }
    
#define __UNIFRAC 1
#endif
