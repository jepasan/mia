#include <iostream>
#include <Rcpp.h>
#include <vector>
#include "api_s.hpp"
#include "tree_s.hpp"

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
Rcpp::List faith_pd_new(const Rcpp::S4 & treeSE, Rcpp::String tree){
    r_vec* result = NULL;
    ComputeStatus status;
    std::string newick(tree.get_cstring());
    status = faith_pd_one_off(treeSE, &result, tree);
    vector<double> values;
    for(int i = 0; i < result->n_samples; i++){
        values.push_back(result->values[i]);
    }
    
    Rcpp::List rlist = Rcpp::List::create(Rcpp::Named("n_samples") = result->n_samples,
                                          Rcpp::Named("faith_pd") = values);
    
    destroy_results_vec(result);
     
     return rlist;
    
    /*
    Rcpp::List rowTree = treeSE.slot("rowTree");
    const List & phylo = rowTree["phylo"];
    const Rcpp::NumericMatrix & edge = phylo["edge"];
    
    return rowTree;
    */
}


std::vector<bool> newick(Rcpp::String ins) {
    std::string newick(ins.get_cstring());
    std::vector<bool> result = std::vector<bool>();
    char last_structure;
    bool potential_single_descendent = false;
    int count = 0;
    bool in_quote = false;
    for(auto c = newick.begin(); c != newick.end(); c++) {
        if(*c == '\'') 
            in_quote = !in_quote;
        
        if(in_quote)
            continue;
        
        switch(*c) {
        case '(':
            // opening of a node
            count++;
            result.push_back(true);
            last_structure = *c;
            potential_single_descendent = true;
            break;
        case ')':
            // closing of a node
            if(potential_single_descendent || (last_structure == ',')) {
                // we have a single descendent or a last child (i.e. ",)" scenario)
                count += 3;
                result.push_back(true);
                result.push_back(false);
                result.push_back(false);
                potential_single_descendent = false;
            } else {
                // it is possible still to have a single descendent in the case of 
                // multiple single descendents (e.g., (...()...) )
                count += 1;
                result.push_back(false);
            }
            last_structure = *c;
            break;
        case ',':
            if(last_structure != ')') {
                // we have a new tip
                count += 2;
                result.push_back(true);
                result.push_back(false);
            }
            potential_single_descendent = false;
            last_structure = *c;
            break;
        default:
            break;
        }
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::LogicalVector treetest(std::string n){
    
    Rcpp::LogicalVector result = Rcpp::LogicalVector();
    std::vector<bool> raw = newick(n);
    if(raw.size() > 0) {
        std::cout << raw.size();
    }
    result = Rcpp::LogicalVector::import(raw.begin(), raw.end());
    return result;
    
}



