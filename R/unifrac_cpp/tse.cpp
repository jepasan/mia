/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <vector>

#include "tse.hpp"

#include <Rcpp.h>

using namespace su;

tse::tse(const Rcpp::S4 & treeSE) {
    sample_ids = std::vector<std::string>(); 
    obs_ids = std::vector<std::string>();
    
    Rcpp::S4 colData = treeSE.slot("colData");
    Rcpp::StringVector rownames = colData.slot("rownames");
    sample_ids = Rcpp::as<std::vector<std::string>>(rownames);
    
    Rcpp::List rowTree = treeSE.slot("rowTree");
    Rcpp::List phylo = rowTree["phylo"];
    Rcpp::StringVector tip_label = phylo["tip.label"];
    obs_ids = Rcpp::as<std::vector<std::string>>(tip_label);
    
    Rcpp::S4 assays = treeSE.slot("assays");
    Rcpp::S4 data = assays.slot("data");
    Rcpp::List listData = data.slot("listData");
    assay = Rcpp::as<Rcpp::NumericMatrix>(listData["counts"]);

    n_samples = sample_ids.size();
    n_obs = obs_ids.size();

    /* define a mapping between an ID and its corresponding offset */
    obs_id_index = std::unordered_map<std::string, uint32_t>();
    sample_id_index = std::unordered_map<std::string, uint32_t>();

    create_id_index(obs_ids, obs_id_index);
    create_id_index(sample_ids, sample_id_index);
    
    sample_counts = get_sample_counts();
    
}

tse::~tse() {
    
}

void tse::create_id_index(std::vector<std::string> &ids, 
                           std::unordered_map<std::string, uint32_t> &map) {
    uint32_t count = 0;
    map.reserve(ids.size());
    for(auto i = ids.begin(); i != ids.end(); i++, count++) {
        map[*i] = count;
    }
}


std::vector<double> tse::get_obs_data(const std::string &id) const {
    std::vector<double> out = std::vector<double>();
    uint32_t idx = obs_id_index.at(id);
    for(unsigned int i = 0; i < n_samples; i++) {
        out.push_back(assay(idx, i));
    }
    return out;
}

std::vector<double> tse::get_sample_counts() {
    std::vector<double> sample_counts = std::vector<double>();
    
    for(unsigned int i = 0; i < n_samples; i++) {
        unsigned int sum = 0;
        for(unsigned int j = 0; j < n_obs; j++){
            sum += assay(j, i);
        }
        sample_counts.push_back(sum);
    }
    return(sample_counts);
}
