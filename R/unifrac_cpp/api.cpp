#include "api.hpp"
#include "biom.hpp"
#include "tree.hpp"
#include "unifrac.hpp"

#include <fstream>
#include <iomanip>
#include <thread>
#include <cstring>
#include <stdlib.h> 

#include <fcntl.h>
#include <unistd.h>

#include <Rcpp.h>

using namespace su;
using namespace std;

std::vector<double> faith_pd_one_off(const Rcpp::S4 & treeSE, bool isRooted){
    
    // Check that tree and table are non-empty and match before calling the c++ code
    // shear the tree (to contain only the obs in the table?) - Also should be done before the call?
    
    std::cout << "Start\n";
    su::BPTree tree = su::BPTree(treeSE, isRooted);      
    std::cout << "Tree ok\n";
    su::tse table = su::tse(treeSE);
    std::cout << "Table ok\n";

    std::vector<double> results =  su::faith_pd(table, tree);
    std::cout << "Results ok\n";

    // compute faithpd
    return results;
    //return std::vector<double>();
}