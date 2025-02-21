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
    su::Assay table = su::Assay(assay);
    
    std::unordered_set<std::string> to_keep(table.obs_ids.begin(),table.obs_ids.end());
    su::BPTree tree_sheared = tree.shear(to_keep).collapse();
    
    su::PropStack propstack(table.n_samples);
    
    uint32_t node;
    std::vector<double> node_proportions;
    double length;
    
    std::vector<double> results = std::vector<double>(table.n_samples, 0.0);
    
    // for node in postorderselect
    const unsigned int max_k = (tree_sheared.nparens>1) ? ((tree_sheared.nparens / 2) - 1) : 0;
    for(unsigned int k = 0; k < max_k; k++) {
        node = tree_sheared.postorderselect(k);
        
        // get branch length
        length = tree_sheared.lengths[node];
        
        // get node proportions and set intermediate scores
        node_proportions = set_proportions(tree_sheared, node, table, propstack, false); // this would probably be the most likely culprit for something going wrong
        
        for (unsigned int sample = 0; sample < table.n_samples; sample++){
            // calculate contribution of node to score
            // Is it possible to somehow set the proportions to 0 if we're dealing with the root in a include.root=FALSE scenario?
            results[sample] += (node_proportions[sample] > 0) * length;
        }
    }
    
    Rcpp::NumericVector faith = Rcpp::NumericVector(results.size());
    
    for(unsigned int i = 0; i < results.size(); i++){
        faith[i] = results[i];
    }
    
    return faith;
}

// [[Rcpp::export]]
double sumrowt(const Rcpp::List & rowTree){
  
  Rcpp::NumericVector edgelength = rowTree["edge.length"];
  
  const uint32_t n_edges = edgelength.size();
  
  //Used to find the correct lengths for the nodes - Includes the root
  double edge_v = 0.0;
  
  for(unsigned int i = 0; i < n_edges; i++){
    edge_v += edgelength[i];
  }
  
  return edge_v;
}