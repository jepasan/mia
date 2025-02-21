/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "tree.h"
#include "assay.h"
#include "propmap.h"

#include <cstdlib>
#include <thread>
#include <signal.h>
#include <stdarg.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>

#include <Rcpp.h>

using namespace su;

PropMap::PropMap(uint32_t vecsize) 
: prop_map()
, defaultsize(vecsize)
{
    prop_map.reserve(1000);
}

PropMap::~PropMap() {
}

std::vector<double> PropMap::get(uint32_t i) {
    if(prop_map.count(i) > 0){
        return prop_map.at(i);
    }
    else {
        return(std::vector<double>());  
    } 
}

void PropMap::clear(uint32_t i) {
    prop_map[i] = std::vector<double>();
}

void PropMap::update(uint32_t node, std::vector<double> vec) {
    prop_map[node] = vec;
}

std::vector<double> su::set_proportions(const BPTree &tree,
                         uint32_t node,
                         const Assay &table,
                         PropMap &ps,
                         bool normalize) {
    
    std::vector<double> props = std::vector<double>();
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
        
        while(current <= right && current != 0) {
            std::vector<double> vec = ps.get(current);  // pull from prop map
            ps.clear(current);  // remove from prop map, place back on stack
            
            for(unsigned int i = 0; i < table.n_samples; i++)
                props[i] = props[i] + vec[i];
            
            current = tree.rightsibling(current);
        }
        
        //std::cout << "n " << props[0] << " " << props[1] << " " << props[2] << "\n";
        
    }
    ps.update(node, props);
    return(props);
}
