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
#include <cstdlib>
#include <thread>
#include <signal.h>
#include <stdarg.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>

#include "unifrac_internal_s.hpp"

using namespace su;

template<class TFloat>
PropStack<TFloat>::PropStack(uint32_t vecsize) 
: prop_stack()
, prop_map()
, defaultsize(vecsize)
{
    prop_map.reserve(1000);
}

template<class TFloat>
PropStack<TFloat>::~PropStack() {
    // drain stack
    for(unsigned int i = 0; i < prop_stack.size(); i++) {
        TFloat *vec = prop_stack.top();
        prop_stack.pop();
        free(vec);
    }

    // drain the map
    for(auto it = prop_map.begin(); it != prop_map.end(); it++) {
        TFloat *vec = it->second;
        free(vec);
    }
    prop_map.clear();
}

template<class TFloat>
TFloat* PropStack<TFloat>::get(uint32_t i) {
    return prop_map[i];
}

template<class TFloat>
void PropStack<TFloat>::push(uint32_t node) {
    TFloat* vec = prop_map[node];
    prop_map.erase(node);
    prop_stack.push(vec);
}

template<class TFloat>
TFloat* PropStack<TFloat>::pop(uint32_t node) {
    /*
     * if we don't have any available vectors, create one
     * add it to our record of known vectors so we can track our mallocs
     */
    void *vec;
    int err = 0;
    if(prop_stack.empty()) {
        // Linux-specific code
        /*err = posix_memalign((void **)&vec, 32, sizeof(TFloat) * defaultsize);
        if(vec == NULL || err != 0) {
            fprintf(stderr, "Failed to allocate %zd bytes, err %d; [%s]:%d\n",
                    sizeof(TFloat) * defaultsize, err, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }*/
        // Windows-specific code
        vec = _aligned_malloc(sizeof(TFloat) * defaultsize, 32);
        if(vec == NULL || !vec) {
            _get_errno(&err);
            fprintf(stderr, "Failed to allocate %zd bytes, err %d; [%s]:%d\n",
                    sizeof(TFloat) * defaultsize, err, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }
    else {
        vec = prop_stack.top();
        prop_stack.pop();
    }

    prop_map[node] = (TFloat*) vec;
    return (TFloat*) vec;
}

// make sure they get instantiated
template class su::PropStack<float>;
template class su::PropStack<double>;


template<class TFloat>
void su::set_proportions(TFloat* __restrict__ props,
                         const BPTree &tree,
                         uint32_t node,
                         const biom_interface &table,
                         PropStack<TFloat> &ps,
                         bool normalize) {
    if(tree.isleaf(node)) {
       table.get_obs_data(tree.names[node], props); // Here we basically just need the row for the specified node
       if (normalize) {
#pragma omp parallel for schedule(static)
        for(unsigned int i = 0; i < table.n_samples; i++) {
           props[i] /= table.sample_counts[i];
        }
       }

    } else {
        unsigned int current = tree.leftchild(node);
        unsigned int right = tree.rightchild(node);

#pragma omp parallel for schedule(static)
        for(unsigned int i = 0; i < table.n_samples; i++)
            props[i] = 0;

        while(current <= right && current != 0) {
            TFloat * __restrict__ vec = ps.get(current);  // pull from prop map
            ps.push(current);  // remove from prop map, place back on stack

#pragma omp parallel for schedule(static)
            for(unsigned int i = 0; i < table.n_samples; i++)
                props[i] = props[i] + vec[i];

            current = tree.rightsibling(current);
        }
    }
}

// make sure they get instantiated
template void su::set_proportions(float* __restrict__ props,
                                  const BPTree &tree,
                                  uint32_t node,
                                  const biom_interface &table,
                                  PropStack<float> &ps,
                                  bool normalize);
template void su::set_proportions(double* __restrict__ props,
                                  const BPTree &tree,
                                  uint32_t node,
                                  const biom_interface &table,
                                  PropStack<double> &ps,
                                  bool normalize);
