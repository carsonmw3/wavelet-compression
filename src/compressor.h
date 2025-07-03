#pragma once

#include "box-structs.h"

// Compresses a multibox, and writes the result to a file. First,
// performs wavelet decomposition on the box, finds a coefficient
// threshold based on "keep", rle-encodes the remainin coefficients,
// serializes them, and writes them to files using lzma.
std::vector<CompressedWavelet> compress(multiBox3D&      box,
                           std::vector<int> components,
                           double           keep,
                           int              time,
                           int              level,
                           int              box_index,
                           std::string      compressed_dir);

