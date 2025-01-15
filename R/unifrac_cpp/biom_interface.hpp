/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2021-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */


#ifndef _UNIFRAC_BIOM_INTERFACE_H
#define _UNIFRAC_BIOM_INTERFACE_H

#include <vector>
#include <string>

#include <Rcpp.h>

//Faith calculations mainly need n_samples, get_obs_data and sample_counts
//sample_counts - OK
//n_samples - OK
//get_obs_data

namespace su {
    class tse_interface {
        public:
            // cache the IDs contained within the table
            std::vector<std::string> sample_ids;
            std::vector<std::string> obs_ids;

            uint32_t n_samples;  // the number of samples
            uint32_t n_obs;      // the number of observations
            std::vector<double> sample_counts; // Counts summed per sample

            /* default constructor
             *
             * Automatically create the needed objects.
             * All other initialization happens in children constructors.
             */
            tse_interface() {}

            /* default destructor
             *
             * Automatically destroy the objects.
             * All other cleanup must have been performed by the children constructors.
             */
            virtual ~tse_interface() {}

            /* get a dense vector of observation data
             *
             * @param id The observation ID to fetch
             * @param out An allocated array of at least size n_samples. 
             *      Values of an index position [0, n_samples) which do not
             *      have data will be zero'd.
             */
            virtual std::vector<double> get_obs_data(const std::string &id) const = 0;
      };
}

#endif /* _UNIFRAC_BIOOM_INTERFACE_H */
