#include <iostream>
#include <vector>

#include <Rcpp.h>

#include "assay.hpp"
#include "tree.hpp"
#include "propstack.hpp"


/* Access the C++ implementation of the fast Faith's PD algorithm from R
 *  
 * treeSE <const Rcpp::S4 &> an R TreeSummarizedExperiment object.
 * faith <Rcpp::NumericVector> the resulting vector of computed Faith PD values
 *
 * This functions makes several assumptions about treeSE:
 * 
 *  - It must contain a non-empty counts assay and a non-empty RowTree
 *  - The RowTree must be sorted in cladewise order
 *  - The RowTree must be rooted - Unrooted trees can be passed without error, but don't produce correct results
 *  
 // Check that tree and table are non-empty and match before calling the c++ code
 // shear the tree (to contain only the obs in the table?) - Also should be done before the call?
 // Assure that tree does not contain ids that are not in table
 * 
 */
// [[Rcpp::export]]
Rcpp::NumericVector faith_cpp(const Rcpp::NumericMatrix & assay, const Rcpp::List & rowTree){
    
    su::BPTree tree = su::BPTree(rowTree);      
    su::Assay table = su::Assay(assay, rowTree);
    
    su::PropStack propstack(table.n_samples);
    
    uint32_t node;
    std::vector<double> node_proportions;
    double length;
    
    std::vector<double> results = std::vector<double>(table.n_samples, 0.0);
    
    // for node in postorderselect
    for(unsigned int k = 0; k < (tree.nparens / 2) - 1; k++) {
        node = tree.postorderselect(k);
        
        // get branch length
        length = tree.lengths[node];
        
        // get node proportions and set intermediate scores
        node_proportions = set_proportions(tree, node, table, propstack); // this would probably be the most likely culprit for something going wrong
        
        for (unsigned int sample = 0; sample < table.n_samples; sample++){
            // calculate contribution of node to score
            // Is it possible to somehow set the proportions to 0 if we're dealing with the root in a include.root=FALSE scenario?
            //if(sample == 0) std::cout << k << " " << node_proportions[sample] << " " << (node_proportions[sample] > 0) << " " << length << "\n";
            results[sample] += (node_proportions[sample] > 0) * length;
        }
    }
    
    Rcpp::NumericVector faith = Rcpp::NumericVector(results.size());
    
    for(unsigned int i = 0; i < results.size(); i++){
        faith[i] = results[i];
    }
    
    return faith;
}