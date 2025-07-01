#pragma once

#include "box-structs.h"

// Decompresses one box at the given file path
Box3D decompress (std::string file_path,
                  int         time,
                  int         level,
                  int         component,
                  int         box_idx);

// Takes a compressed box written in raw bytes and deserializes it back into
// a CompressedWavelet structure for decoding/decomposition
CompressedWavelet deserialize_compressed_wavelet(const std::string& data);

// "Undoes"/decomposes a flattened vector of coefficients from a Haar wavelet
// decomposition for a box, and returns the decomposed box
Box3D inverse_wavelet_decompose(std::vector<float> flat, int x, int y, int z);
