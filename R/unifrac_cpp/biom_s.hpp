/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */


#ifndef _UNIFRAC_BIOM_H
#define _UNIFRAC_BIOM_H

#include <H5Cpp.h>
#include <H5Dpublic.h>
#include <vector>
#include <unordered_map>

#include "biom_interface_s.hpp"

#include <Rcpp.h>

namespace su {
    class tse : public tse_interface {
        public:
            /* default constructor
             *
             * @param treeSE An R TreeSummarizedExperiment object
             */
            tse(const Rcpp::S4 & treeSE);

            /* default destructor
             *
             * Temporary arrays are freed
             */
            virtual ~tse();

            /* get a dense vector of observation data
             *
             * @param id The observation ID to fetch
             * @param out An allocated array of at least size n_samples. 
             *      Values of an index position [0, n_samples) which do not
             *      have data will be zero'd.
             */
            std::vector<double> get_obs_data(const std::string &id) const;

        private:
            Rcpp::NumericMatrix assay; // Access to the raw sample counts in R's memory

            std::vector<double> get_sample_counts();

            /* At construction, lookups mapping IDs -> index position within an
             * axis are defined
             */
            std::unordered_map<std::string, uint32_t> obs_id_index;
            std::unordered_map<std::string, uint32_t> sample_id_index;
 
            /* load ids from an axis
             *
             * @param path The dataset path to the ID dataset to load
             * @param ids The variable representing the IDs to load into
             */          
            void load_ids(const char *path, std::vector<std::string> &ids);

            /* load the index pointer for an axis
             *
             * @param path The dataset path to the index pointer to load
             * @param indptr The vector to load the data into
             */
            void load_indptr(const char *path, std::vector<uint32_t> &indptr);

            /* count the number of nonzero values and set nnz */
            void set_nnz();

            /* create an index mapping an ID to its corresponding index 
             * position.
             *
             * @param ids A vector of IDs to index
             * @param map A hash table to populate
             */
            void create_id_index(std::vector<std::string> &ids, 
                                 std::unordered_map<std::string, uint32_t> &map);

            // templatized version
            template<class TFloat> std::vector<TFloat> get_obs_data_TT(const std::string &id, TFloat t) const;
            
     };
}

#endif /* _UNIFRAC_BIOM_H */

