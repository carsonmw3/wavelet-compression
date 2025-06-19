#pragma once

#include "box-structs.h"

// TODO: If you are repeating parameters or code more than two-ish times,
// that typically means you'll want to extract them into a struct

void write_loc_dim_to_bin(LocDimData  data,
                          std::string path,
                          std::string out_file,
                          int         num_times,
                          int         num_levels);


LocDimData read_loc_dim_from_bin(std::string const&            path,
                                 std::string const&            in_file,
                                 std::vector<std::vector<int>> counts,
                                 int                           num_times,
                                 int                           num_levels);


void write_box_counts(std::vector<std::vector<int>> counts,
                      std::string const&            path,
                      std::string const&            out_file,
                      int                           num_times,
                      int                           num_levels);


std::vector<std::vector<int>> read_box_counts(std::string path,
                                              std::string in_file,
                                              int         num_times,
                                              int         num_levels);

void write_box_with_loc_dim(Volume3D           volume,
                            std::string        path,
                            std::string        out_file,
                            std::vector<float> location,
                            std::vector<float> dimension);

