#pragma once

#include "box-structs.h"

std::vector<CompressedWavelet> compress(multiBox3D&      box,
                           std::vector<int> components,
                           double           keep,
                           int              time,
                           int              level,
                           int              box_index,
                           std::string      compressed_dir);

