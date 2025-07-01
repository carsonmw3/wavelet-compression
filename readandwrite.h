#pragma once

#include "box-structs.h"
#include "iterator.h"

// writes location or dimension data for a compression run to a binary file
void write_loc_dim_to_bin(LocDimData  data,
                          std::string path,
                          std::string out_file,
                          AMRIterator iterator);


// reads location or dimension data for a compression run from a binary file
LocDimData read_loc_dim_from_bin(std::string const&            path,
                                 std::string const&            in_file,
                                 std::vector<std::vector<int>> counts,
                                 AMRIterator                   iterator,
                                 int                           num_times,
                                 int                           num_levels);


// writes box count data for a compression run to a binary file
void write_box_counts(std::vector<std::vector<int>> counts,
                      std::string const&            path,
                      std::string const&            out_file,
                      int                           num_times,
                      int                           num_levels);


// reads box count data for a compression run from a binary file
std::vector<std::vector<int>> read_box_counts(std::string path,
                                              std::string in_file,
                                              int         num_times,
                                              int         num_levels);

// writes amrex formatting info for a compression run to a binary file
void write_amrexinfo(AMReXInfo   info,
                     std::string path,
                     std::string out_file);

// reads amrex formatting info for a compression run from a binary file
AMReXInfo read_amrex_info(std::string path,
                          std::string in_file);
