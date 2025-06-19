#pragma once

#include "box-structs.h"

#include <utility>

CompressedWavelet read_compressed_wavelet(const std::string& filename);

std::vector<float> rle_decode(std::vector<std::pair<int, float>> rle_encoded,
                              int total_length);

Box3D inverse_wavelet_decompose(std::vector<float> flat, int x, int y, int z);

