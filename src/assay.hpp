/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef __FAITH_ASSAY_H
#define __FAITH_ASSAY_H 1

#include <vector>
#include <unordered_map>

#include <Rcpp.h>

namespace su {
    class Assay {
        public:
            // cache the IDs contained within the table
            std::vector<std::string> sample_ids;
            std::vector<std::string> obs_ids;
            
            uint32_t n_samples;  // the number of samples
            uint32_t n_obs;      // the number of observations
            std::vector<double> sample_counts; // Counts summed per sample
            
            /* default constructor
             *
             * @param treeSE An R TreeSummarizedExperiment object
             */
            Assay(const Rcpp::NumericMatrix & assay, const Rcpp::List & rowTree);

            /* default destructor
             *
             * Temporary arrays are freed
             */
            ~Assay();

            /* get a dense vector of observation data
             *
             * @param id The observation ID to fetch
             * @param out An allocated array of at least size n_samples. 
             *      Values of an index position [0, n_samples) which do not
             *      have data will be zero'd.
             */
            std::vector<double> get_obs_data(const std::string &id) const;

        private:
            Rcpp::NumericMatrix table; // Access to the raw sample counts in R's memory

            std::vector<double> get_sample_counts();

            /* At construction, lookups mapping IDs -> index position within an
             * axis are defined
             */
            std::unordered_map<std::string, uint32_t> obs_id_index;
            std::unordered_map<std::string, uint32_t> sample_id_index;

            /* create an index mapping an ID to its corresponding index 
             * position.
             *
             * @param ids A vector of IDs to index
             * @param map A hash table to populate
             */
            void create_id_index(std::vector<std::string> &ids, 
                                 std::unordered_map<std::string, uint32_t> &map);
     };
}

#endif /* __FAITH_ASSAY_H */

