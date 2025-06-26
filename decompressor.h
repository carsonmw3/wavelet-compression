#pragma once

#include "box-structs.h"

Box3D decompress (std::string file_path,
                  int         time,
                  int         level,
                  int         component,
                  int         box_idx);

CompressedWavelet deserialize_compressed_wavelet(const std::string& data);

Box3D inverse_wavelet_decompose(std::vector<float> flat, int x, int y, int z);
