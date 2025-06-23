#pragma once

#include "box-structs.h"

void write_plotfiles(std::vector<std::vector<std::vector<multiBox3D>>> &data,
                     LocDimData                                        locations,
                     LocDimData                                        dimensions,
                     int                                               num_times,
                     int                                               num_levels,
                     int                                               num_components,
                     AMReXInfo                                         amrexinfo,
                     std::string                                       out);

