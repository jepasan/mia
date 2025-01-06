#include <iostream>
#include <Rcpp.h>
#include <vector>
#include "api_s.hpp"
#include "tree_s.hpp"

#include <stack>

using namespace std;
using namespace Rcpp;


/*
// [[Rcpp::export]]
Rcpp::List faith_pd(const char* table, const char* tree){
    r_vec* result = NULL;
    ComputeStatus status;
    status = faith_pd_one_off(table, tree, &result);
    vector<double> values;
    for(int i = 0; i < result->n_samples; i++){
        values.push_back(result->values[i]);
    }
    
    return Rcpp::List::create(Rcpp::Named("n_samples") = result->n_samples,
                              Rcpp::Named("faith_pd") = values);
    
}
*/

// [[Rcpp::export]]
void faith_pd_new(const Rcpp::S4 & treeSE){
    
    std::vector<double> results = faith_pd_one_off(treeSE);
    
    std::cout << results.size() << "\n";
    
    if(results.size() >= 20){
        for(unsigned int i = 0; i < 20; i++){
            std::cout << results[i] << "\n";
        }
    }
    //get_sample_counts(treeSE);
    /*
    r_vec* result = NULL;
    ComputeStatus status;
    status = faith_pd_one_off(treeSE, &result);
    
    vector<double> values;
    for(int i = 0; i < result->n_samples; i++){
        values.push_back(result->values[i]);
    }
    
    Rcpp::List rlist = Rcpp::List::create(Rcpp::Named("n_samples") = result->n_samples,
                                          Rcpp::Named("faith_pd") = values);
    
    destroy_results_vec(&result);
    
    return rlist;
    */
    
    //Rcpp::List rowTree = treeSE.slot("rowTree");
    //const List & phylo = rowTree["phylo"];
    //const Rcpp::NumericMatrix & edge = phylo["edge"];
    
    //return rowTree;
}


