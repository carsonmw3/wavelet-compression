#pragma once

#include "grid.h"
#include <string>

// A 1-component box of data
using Box3D      = Grid3D<float>;

// A multi-component box of data (stored as a vector of 1-component boxes)
using multiBox3D = std::vector<Box3D>;

// Stores the location of a box (where in the volume it's indexed from)
using Location   = std::vector<int>;

// Stores the dimensions of a box
using Dimensions = std::vector<int>;

// location or dimension data for all boxes in a compression run
using LocDimData = std::vector<std::vector<std::vector<std::vector<int>>>>;

// Stores info about the compression run for use in decompression
struct RunInfo {
    int              min_time;
    int              max_time;
    int              min_level;
    int              max_level;
    std::vector<int> components;
};

// Stores relevant data for one level at one timestep
struct LevelData {
    std::vector<multiBox3D> boxes;
    std::vector<Location>   locations;
    std::vector<Dimensions> dimensions;
    int                     box_count;
    std::vector<float>      min_values;
    std::vector<float>      max_values;
};

// Stores relevant data about the AMReX parameters needed when
// writing plotfiles from decompressed data
struct AMReXInfo {
    std::vector<std::string>         comp_names;
    std::vector<std::vector<double>> geomcellinfo;
    std::vector<int>                 ref_ratios;
    std::vector<long double>         true_times;
    std::vector<std::vector<int>>    level_steps;
    int                              xDim;
    int                              yDim;
    int                              zDim;
};

// Stores all relevant data for a single compression run
struct AllData {
    std::vector<std::vector<std::vector<multiBox3D>>> boxes;
    LocDimData                                        locations;
    LocDimData                                        dimensions;
    std::vector<std::vector<int>>                     box_counts;
    std::vector<float>                                min_values;
    std::vector<float>                                max_values;
    AMReXInfo                                         amrexinfo;
    RunInfo                                           runinfo;
};

// Stores relevant compressed data for a single box
struct CompressedWavelet {
    std::vector<int>                     shape;
    std::vector<int>                     coeff_shape;
    std::vector<std::pair<int, float>>   rle_encoded;
    bool                                 need32;
};

