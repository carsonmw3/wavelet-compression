#pragma once

#include "grid.h"
#include <cstdint>
#include <utility>

using Box3D      = Grid3D<float>;
using Volume3D   = Grid3D<float>;
using Location   = std::vector<int>;
using Dimensions = std::vector<int>;
using LocDimData = std::vector<std::vector<std::vector<std::vector<int>>>>;
// TODO: ^ could we use a Grid< Vec<Grid<int>> > or similar here?

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
