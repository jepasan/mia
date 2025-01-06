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
#include <Rcpp.h>

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
}

template<class TFloat>
std::vector<TFloat> PropStack<TFloat>::get(uint32_t i) {
    if(prop_map.count(i) > 0){
        return prop_map.at(i);
    }
    else {
        return(std::vector<TFloat>());  
    } 
}

template<class TFloat>
void PropStack<TFloat>::clear(uint32_t i) {
    prop_map[i] = std::vector<TFloat>();
}

template<class TFloat>
void PropStack<TFloat>::update(uint32_t node, std::vector<TFloat> vec) {
    prop_map[node] = vec;
}

// make sure they get instantiated
template class su::PropStack<float>;
template class su::PropStack<double>;


template<class TFloat>
std::vector<TFloat> su::set_proportions(const BPTree &tree,
                         uint32_t node,
                         const tse_interface &table,
                         PropStack<TFloat> &ps,
                         bool normalize) {
    
    std::vector<TFloat> props = std::vector<TFloat>();
    if(tree.isleaf(node)) {
        std::string leaf = tree.names[node];
        props = table.get_obs_data(leaf); // Here we basically just need the row for the specified node
        if (normalize) {
//#pragma omp parallel for schedule(static)
            for(unsigned int i = 0; i < table.n_samples; i++) {
               props[i] /= table.sample_counts[i];
            }
        } 
    } else {
        unsigned int current = tree.leftchild(node);
        unsigned int right = tree.rightchild(node);

//#pragma omp parallel for schedule(static)
        
        for(unsigned int i = 0; i < table.n_samples; i++){
            props.push_back(0);
        }
        
        while(current <= right && current != 0) {
            std::vector<TFloat> vec = ps.get(current);  // pull from prop map
            ps.clear(current);  // remove from prop map, place back on stack
//#pragma omp parallel for schedule(static)
            for(unsigned int i = 0; i < table.n_samples; i++)
                props[i] = props[i] + vec[i];
            
            current = tree.rightsibling(current);
        }
        
    }
    ps.update(node, props);
    return(props);
}

// make sure they get instantiated
template std::vector<double> su::set_proportions(const BPTree &tree,
                                  uint32_t node,
                                  const tse_interface &table,
                                  PropStack<double> &ps,
                                  bool normalize);
