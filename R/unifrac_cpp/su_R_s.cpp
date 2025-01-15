#include <iostream>
#include <Rcpp.h>
#include <vector>
#include "api_s.hpp"
#include "tree_s.hpp"

#include <stack>

// [[Rcpp::export]]
Rcpp::NumericVector faith_pd(const Rcpp::S4 & treeSE, bool isRooted){
    
    std::vector<double> results = faith_pd_one_off(treeSE, isRooted);
    
    Rcpp::NumericVector faith = Rcpp::NumericVector(results.size());
    
    for(unsigned int i = 0; i < results.size(); i++){
        faith[i] = results[i];
    }
    
    return(faith);
}


// [[Rcpp::export]]
Rcpp::LogicalVector rowTree_to_bp(const Rcpp::List & rowTree) {
    Rcpp::NumericMatrix edge = rowTree["edge"];
    Rcpp::StringVector tips = rowTree["tip.label"];
    std::vector<bool> structure = std::vector<bool>();
    
    uint32_t ntips = tips.size(); // phylo tips are always numbered from 1 to number of tips;
    
    std::stack<unsigned int> nodes; // Keeps track of the branch's internal nodes
    
    int currentNode = 0;
    int nextNode = 0;
    
    // Goal: Insert true when a branch starts, a false when it closes, and a true-false for each tip.
    
    for (unsigned int i = 0; i < edge.nrow(); i++){
        currentNode = edge(i, 0);
        nextNode = edge(i, 1);
        
        if(nodes.size() > 0 && currentNode < nodes.top()) {
            // We've exhausted the branch and moved backwards in the tree
            do {
                nodes.pop();
                structure.push_back(false);
            } while(currentNode != nodes.top());
        }
        
        if(nodes.size() == 0 || currentNode > nodes.top() ) {
            // We are either at the root, or entering a new node
            // What if the tree is unrooted?
            nodes.push(currentNode);
            structure.push_back(true);
            
        }
        
        if(nextNode <= ntips) {
            // We've found a tip
            structure.push_back(true);
            structure.push_back(false);
        }
        
        if(i == edge.nrow() - 1) {
            // We've reached the end of the tree
            do {
                nodes.pop();
                structure.push_back(false);
            } while(nodes.size() > 0);
        }
    }
    
    Rcpp::LogicalVector bp = Rcpp::LogicalVector(structure.size());
    
    for(unsigned int i = 0; i < structure.size(); i++){
        bp[i] = structure[i];
    }
    
    return(bp);
}


// [[Rcpp::export]]
Rcpp::LogicalVector newick_to_bp(std::string newick) {
    char last_structure;
    bool potential_single_descendent = false;
    int count = 0;
    bool in_quote = false;
    std::vector<bool> structure;
    for(auto c = newick.begin(); c != newick.end(); c++) {
        if(*c == '\'') 
            in_quote = !in_quote;
        
        if(in_quote)
            continue;
        
        switch(*c) {
        case '(':
            // opening of a node
            count++;
            structure.push_back(true);
            last_structure = *c;
            potential_single_descendent = true;
            break;
        case ')':
            // closing of a node
            if(potential_single_descendent || (last_structure == ',')) {
                // we have a single descendent or a last child (i.e. ",)" scenario)
                count += 3;
                structure.push_back(true);
                structure.push_back(false);
                structure.push_back(false);
                potential_single_descendent = false;
            } else {
                // it is possible still to have a single descendent in the case of 
                // multiple single descendents (e.g., (...()...) )
                count += 1;
                structure.push_back(false);
            }
            last_structure = *c;
            break;
        case ',':
            if(last_structure != ')') {
                // we have a new tip
                count += 2;
                structure.push_back(true);
                structure.push_back(false);
            }
            potential_single_descendent = false;
            last_structure = *c;
            break;
        default:
            break;
        }
    }
    
    Rcpp::LogicalVector bp = Rcpp::LogicalVector(structure.size());
    
    for(unsigned int i = 0; i < structure.size(); i++){
        bp[i] = structure[i];
    }
    
    return(bp);
}
