/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "tree.hpp"
#include "assay.hpp"
#include "propstack.hpp"

#include <cstdlib>
#include <thread>
#include <signal.h>
#include <stdarg.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>

#include <Rcpp.h>

using namespace su;

PropStack::PropStack(uint32_t vecsize) 
: prop_map()
, defaultsize(vecsize)
{
    prop_map.reserve(1000);
}

PropStack::~PropStack() {
}

std::vector<double> PropStack::get(uint32_t i) {
    if(prop_map.count(i) > 0){
        return prop_map.at(i);
    }
    else {
        return(std::vector<double>());  
    } 
}

void PropStack::clear(uint32_t i) {
    prop_map[i] = std::vector<double>();
}

void PropStack::update(uint32_t node, std::vector<double> vec) {
    prop_map[node] = vec;
}

std::vector<double> su::set_proportions(const BPTree &tree,
                         uint32_t node,
                         const Assay &table,
                         PropStack &ps,
                         bool normalize) {
    
    std::vector<double> props = std::vector<double>();
    
    //propstack.clear(node); the current node is popped from propstack at every loop, replacing the vector with an empty one that then gets filled...
  
  
    if(tree.isleaf(node)) {
        std::string leaf = tree.names[node];
      
        props = table.get_obs_data(leaf); // Here we basically just need the row for the specified node
        if (normalize) {
            for(unsigned int i = 0; i < table.n_samples; i++) {
               props[i] /= table.sample_counts[i];
            }
        } 
    } else {
        unsigned int current = tree.leftchild(node);
        unsigned int right = tree.rightchild(node);

        for(unsigned int i = 0; i < table.n_samples; i++){
            props.push_back(0);
        }
        ps.update(node, props);
        
        while(current <= right && current != 0) {
            std::vector<double> vec = ps.get(current);  // pull from prop map
            ps.clear(current);  // remove from prop map, place back on stack
            
            for(unsigned int i = 0; i < table.n_samples; i++)
                props[i] = props[i] + vec[i];
            ps.update(node, props);
            
            current = tree.rightsibling(current);
        }
        
    }
    ps.update(node, props);
    return(props);
}