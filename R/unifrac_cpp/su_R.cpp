#include <iostream>
#include <vector>

#include <Rcpp.h>

#include "api.hpp"
#include "tree.hpp"

// [[Rcpp::export]]
Rcpp::NumericVector faith_pd(const Rcpp::S4 & treeSE, bool isRooted){
    
    std::vector<double> results = faith_pd_one_off(treeSE, isRooted);
    
    Rcpp::NumericVector faith = Rcpp::NumericVector(results.size());
    
    for(unsigned int i = 0; i < results.size(); i++){
        faith[i] = results[i];
    }
    
    return(faith);
}
