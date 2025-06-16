#pragma once

#include "grid.h"

using Box3D      = Grid3D<float>;
using Volume3D   = Grid3D<float>;
using Location   = std::vector<int>;
using Dimensions = std::vector<int>;
using LocDimData = std::vector<std::vector<std::vector<std::vector<int>>>>;

struct LevelData {
    std::vector<Box3D>      boxes;
    std::vector<Location>   locations;
    std::vector<Dimensions> dimensions;
    int                     box_count;
    float                   min_value;
    float                   max_value;
};

struct AllData {
        std::vector<std::vector<std::vector<Box3D>>> boxes;
        LocDimData                                   locations;
        LocDimData                                   dimensions;
        std::vector<std::vector<int>>                box_counts;
        float                                        min_value;
        float                                        max_value;
    };


struct CompressedWavelet {
    std::vector<int> shape;
    std::vector<int> coeff_shape;
    std::vector<std::pair<int, int16_t>> rle_encoded;
};

