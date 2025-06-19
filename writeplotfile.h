#pragma once

#include "box-structs.h"

void write_plotfiles(std::vector<std::vector<std::vector<Box3D>>> &data,
                     LocDimData                                   locations,
                     LocDimData                                   dimensions,
                     int                                          num_times,
                     int                                          num_levels,
                     AMReXInfo                                    amrexinfo,
                     std::string                                  out);

