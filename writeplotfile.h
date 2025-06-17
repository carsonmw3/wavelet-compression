#pragma once

#include "box-structs.h"

void write_plotfiles(std::vector<std::vector<std::vector<Box3D>>> &data,
                     LocDimData                                   locations,
                     LocDimData                                   dimensions,
                     int                                          num_times,
                     int                                          num_levels,
                     std::vector<std::vector<double>>             geomcellinfo,
                     std::vector<std::vector<int>>                refratios,
                     std::vector<long double>                     true_times,
                     std::vector<std::vector<int>>                level_steps,
                     int                                          xDim,
                     int                                          yDim,
                     int                                          zDim,
                     std::string                                  out);

